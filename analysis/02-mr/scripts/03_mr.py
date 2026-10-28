"""Two-sample Mendelian randomisation between bone mineral density / osteoporosis
and cardiovascular endpoints, from GWAS Catalog harmonised summary statistics.

Pairs are run in both directions (bone -> vascular and vascular -> bone). Methods:
inverse-variance weighted (fixed and random effects), MR-Egger, weighted median,
weighted mode and simple median, with Cochran's Q, the MR-Egger intercept,
leave-one-out and Steiger filtering as sensitivity analyses.

Usage (project .venv):
    python analysis/02-mr/scripts/03_mr.py --self-test
    python analysis/02-mr/scripts/03_mr.py --pairs eBMD:CAD,eBMD:MI --threads 4
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import _sumstats

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data" / "raw" / "gwas"
TAB = ROOT / "results" / "tables" / "02-mr"
CACHE = ROOT / "data" / "processed" / "mr" / "outcome_snps"

FILES = _sumstats.FILES
BINARY = {"CAD", "STROKE", "HF", "MI", "OP"}
RNG = np.random.default_rng(20261021)


# ---------------------------------------------------------------- data plumbing
def load_instruments(key):
    path = TAB / f"instruments_{key}.csv"
    if not path.exists():
        return None
    d = pd.read_csv(path)
    return d.rename(columns={"beta": "bx", "se": "sx", "p": "px", "eaf": "eaf_x", "n": "n_x"})


def outcome_at_snps(key, snps):
    """Effect of each instrument SNP from the outcome file (cached per trait)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    cpath = CACHE / f"{key}.csv"
    got = pd.DataFrame()
    if cpath.exists():
        got = pd.read_csv(cpath)
        if set(snps).issubset(set(got["SNP"])):
            return got[got["SNP"].isin(snps)]
    # path_of knows that FinnGen files live outside data/raw/gwas
    if not _sumstats.path_of(key).exists():
        return None
    new = _sumstats.read_at_snps(key, snps)
    if new is None or new.empty:
        return None
    if len(got):
        new = pd.concat([got[~got["SNP"].isin(new["SNP"])], new], ignore_index=True)
    new.to_csv(cpath, index=False)
    return new[new["SNP"].isin(snps)]


def harmonise(exp, out):
    """Align the outcome to the exposure effect allele and drop unusable SNPs.

    `exp` is an instrument table (bx, sx, px, eaf_x, n_x, effect_allele, other_allele);
    `out` comes from outcome_at_snps (beta, se, p, eaf, n, effect_allele, other_allele).
    """
    e = exp.rename(columns={"effect_allele": "ea_x", "other_allele": "oa_x"})
    o = out.rename(columns={"effect_allele": "ea_y", "other_allele": "oa_y", "eaf": "eaf_y",
                            "beta": "by", "se": "sy", "p": "py", "n": "n_y"})
    m = e.merge(o[["SNP", "ea_y", "oa_y", "eaf_y", "by", "sy", "py", "n_y"]], on="SNP")
    keep = []
    for _, r in m.iterrows():
        a1, a2 = r["ea_x"], r["oa_x"]
        b1, b2 = r["ea_y"], r["oa_y"]
        by, eaf_y = r["by"], r["eaf_y"]
        if (a1, a2) == (b1, b2):
            pass
        elif (a1, a2) == (b2, b1):
            by, eaf_y = -by, None if pd.isna(eaf_y) else 1 - eaf_y
        else:
            continue
        if {a1, a2} in ({"A", "T"}, {"C", "G"}) and not pd.isna(eaf_y) and 0.42 < eaf_y < 0.58:
            continue
        keep.append((r["SNP"], a1, a2, r["eaf_x"], r["bx"], r["sx"], r["px"], r["n_x"],
                     by, r["sy"], r["py"], r["n_y"]))
    return pd.DataFrame(keep, columns=["SNP", "effect_allele", "other_allele", "eaf", "bx", "sx", "px",
                                       "n_x", "by", "sy", "py", "n_y"])


# ---------------------------------------------------------------- MR estimators
def wald(bx, by, sx, sy):
    ratio = by / bx
    se = sy / np.abs(bx)
    return ratio, se


def ivw(bx, by, sx, sy, random=True):
    ratio, se = wald(bx, by, sx, sy)
    w = 1 / se ** 2
    b = np.sum(w * ratio) / np.sum(w)
    se_fixed = np.sqrt(1 / np.sum(w))
    q = np.sum(w * (ratio - b) ** 2)
    dfree = len(ratio) - 1
    p_q = 1 - stats.chi2.cdf(q, dfree) if dfree > 0 else np.nan
    if random and dfree > 0 and q > dfree:
        se_use = se_fixed * np.sqrt(q / dfree)
    else:
        se_use = se_fixed
    p = 2 * stats.norm.sf(abs(b / se_use))
    return dict(b=b, se=se_use, p=p, q=q, q_df=dfree, q_p=p_q, se_fixed=se_fixed)


def egger(bx, by, sx, sy):
    ratio, se = wald(bx, by, sx, sy)
    w = 1 / se ** 2
    X = np.column_stack([np.ones_like(bx), bx])
    W = np.diag(w)
    XtWX = X.T @ W @ X
    beta = np.linalg.solve(XtWX, X.T @ W @ by)
    resid = by - X @ beta
    dfree = len(bx) - 2
    sigma2 = np.sum(w * resid ** 2) / dfree if dfree > 0 else 1.0
    scale = max(1.0, np.sqrt(sigma2))
    cov = np.linalg.inv(XtWX) * sigma2
    se_int, se_slope = np.sqrt(np.diag(cov))
    p_slope = 2 * stats.t.sf(abs(beta[1] / se_slope), dfree) if dfree > 0 else np.nan
    p_int = 2 * stats.t.sf(abs(beta[0] / se_int), dfree) if dfree > 0 else np.nan
    return dict(b=beta[1], se=se_slope, p=p_slope, intercept=beta[0], intercept_se=se_int,
                intercept_p=p_int, scale=scale)


def weighted_median(bx, by, sx, sy, nboot=2000):
    ratio, se = wald(bx, by, sx, sy)
    w = (1 / se ** 2) / np.sum(1 / se ** 2)
    order = np.argsort(ratio)
    r, ww = ratio[order], w[order]
    cum = np.cumsum(ww)
    b = float(r[np.searchsorted(cum, 0.5)])
    boot = []
    k = len(r)
    p = ww
    for _ in range(nboot):
        idx = RNG.choice(k, size=k, replace=True, p=p)
        rb, wb = r[idx], ww[idx]
        o = np.argsort(rb)
        c = np.cumsum(wb[o]) / np.sum(wb)
        boot.append(rb[o][np.searchsorted(c, 0.5)])
    se = float(np.std(boot, ddof=1))
    p_val = 2 * stats.norm.sf(abs(b / se)) if se > 0 else np.nan
    return dict(b=b, se=se, p=p_val)


def weighted_mode(bx, by, sx, sy, phi=1.0, nboot=2000):
    ratio, se = wald(bx, by, sx, sy)
    w = (1 / se ** 2) / np.sum(1 / se ** 2)
    bw = phi * float(np.median(se))
    if bw <= 0:
        bw = 1e-8

    def mode(r, ww):
        def scan(lo, hi, npts):
            grid = np.linspace(lo, hi, npts)
            dens = np.exp(-0.5 * ((grid[:, None] - r[None, :]) / bw) ** 2) @ ww
            return grid, dens
        grid, dens = scan(np.min(r) - 3 * bw, np.max(r) + 3 * bw, 512)
        i = int(np.argmax(dens))
        a, b = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
        grid, dens = scan(a, b, 201)
        return float(grid[int(np.argmax(dens))])

    b = mode(ratio, w)
    boot = []
    k = len(ratio)
    for _ in range(nboot):
        idx = RNG.choice(k, size=k, replace=True, p=w)
        boot.append(mode(ratio[idx], w[idx]))
    se = float(np.std(boot, ddof=1))
    p_val = 2 * stats.norm.sf(abs(b / se)) if se > 0 else np.nan
    return dict(b=b, se=se, p=p_val)


def simple_median(bx, by, sx, sy):
    ratio, _ = wald(bx, by, sx, sy)
    b = float(np.median(ratio))
    se = 1.253 * np.std(ratio, ddof=1) / np.sqrt(len(ratio))
    return dict(b=b, se=float(se), p=2 * stats.norm.sf(abs(b / se)) if se > 0 else np.nan)


def steiger(bx, by, sx, sy, nx, ny):
    """Per-SNP proportion of variance explained, exposure vs outcome (Hemani 2017)."""
    def r2(b, se, n):
        t = b / se
        r = t / np.sqrt(t ** 2 + np.maximum(n - 2, 1))
        return r ** 2, r
    r2x, rx = r2(bx, sx, nx)
    r2y, ry = r2(by, sy, ny)
    sex = 2 * rx * np.sqrt((1 - rx ** 2) ** 2 / np.maximum(nx - 1, 1))
    sey = 2 * ry * np.sqrt((1 - ry ** 2) ** 2 / np.maximum(ny - 1, 1))
    z = (r2x - r2y) / np.sqrt(sex ** 2 + sey ** 2)
    return dict(r2_exposure=float(np.sum(r2x)), r2_outcome=float(np.sum(r2y)),
                steiger_p=float(2 * stats.norm.sf(abs(np.sum(z) / np.sqrt(len(z))))),
                n_wrong_direction=int(np.sum(r2y > r2x)))


def run_pair(exp_key, out_key, min_snps=3, verbose=True):
    inst = load_instruments(exp_key)
    if inst is None or inst.empty:
        return None
    out = outcome_at_snps(out_key, inst["SNP"].tolist())
    if out is None or out.empty:
        return None
    d = harmonise(inst, out)
    if len(d) < min_snps:
        if verbose:
            print(f"  {exp_key} -> {out_key}: only {len(d)} usable SNPs, skipped")
        return None
    bx, by = d["bx"].to_numpy(), d["by"].to_numpy()
    sx, sy = d["sx"].to_numpy(), d["sy"].to_numpy()
    res = []
    for name, fn in [("IVW (fixed)", lambda: {**ivw(bx, by, sx, sy, random=False)}),
                     ("IVW (random)", lambda: ivw(bx, by, sx, sy, random=True)),
                     ("MR-Egger", lambda: egger(bx, by, sx, sy)),
                     ("Weighted median", lambda: weighted_median(bx, by, sx, sy)),
                     ("Weighted mode", lambda: weighted_mode(bx, by, sx, sy)),
                     ("Simple median", lambda: simple_median(bx, by, sx, sy))]:
        r = fn()
        r.update(exposure=exp_key, outcome=out_key, method=name, nsnp=len(d))
        res.append(r)
    st = steiger(bx, by, sx, sy, d["n_x"].to_numpy(), d["n_y"].to_numpy())
    loo = leave_one_out(bx, by, sx, sy)
    summary = dict(exposure=exp_key, outcome=out_key, nsnp=len(d),
                   mean_F=float(d["bx"].abs().div(d["sx"]).pow(2).mean()),
                   min_F=float((d["bx"] / d["sx"]).pow(2).min()),
                   r2_exposure=st["r2_exposure"], r2_outcome=st["r2_outcome"],
                   steiger_p=st["steiger_p"], n_wrong_direction=st["n_wrong_direction"],
                   loo_b_min=loo["b_min"], loo_b_max=loo["b_max"],
                   loo_max_influence=loo["max_shift"])
    if verbose:
        ivw_r = res[1]
        print(f"  {exp_key} -> {out_key}: n={len(d)}  IVW b={ivw_r['b']:.3f} "
              f"({ivw_r['b'] - 1.96 * ivw_r['se']:.3f}, {ivw_r['b'] + 1.96 * ivw_r['se']:.3f}) "
              f"p={ivw_r['p']:.2g}  Q p={ivw_r['q_p']:.2g}  Egger int p={res[2]['intercept_p']:.2g}")
    return pd.DataFrame(res), summary


def leave_one_out(bx, by, sx, sy):
    k = len(bx)
    if k < 4:
        return dict(b_min=np.nan, b_max=np.nan, max_shift=np.nan)
    bs = []
    for i in range(k):
        m = np.ones(k, bool)
        m[i] = False
        bs.append(ivw(bx[m], by[m], sx[m], sy[m], random=True)["b"])
    bs = np.array(bs)
    full = ivw(bx, by, sx, sy, random=True)["b"]
    return dict(b_min=float(bs.min()), b_max=float(bs.max()), max_shift=float(np.max(np.abs(bs - full))))


# ---------------------------------------------------------------- self test
def self_test(n=40, beta_true=0.35, pleiotropy=0.0, seed=7):
    """Simulate a known causal effect and check that the estimators recover it."""
    rng = np.random.default_rng(seed)
    bx = rng.uniform(0.02, 0.12, n) * rng.choice([-1, 1], n)
    sx = rng.uniform(0.005, 0.02, n)
    sy = rng.uniform(0.01, 0.03, n)
    by = beta_true * bx + pleiotropy + rng.normal(0, sy)
    print(f"simulated: true beta={beta_true}, pleiotropy={pleiotropy}, {n} SNPs")
    out = []
    for name, fn in [("IVW (random)", lambda: ivw(bx, by, sx, sy)),
                     ("MR-Egger", lambda: egger(bx, by, sx, sy)),
                     ("Weighted median", lambda: weighted_median(bx, by, sx, sy)),
                     ("Weighted mode", lambda: weighted_mode(bx, by, sx, sy))]:
        r = fn()
        out.append(dict(method=name, b=round(r["b"], 4), se=round(r["se"], 4),
                        p=f"{r['p']:.2g}",
                        extra=f"intercept={r['intercept']:.4f} p={r['intercept_p']:.2g}"
                        if "intercept" in r else ""))
        print(f"  {name:16s} b={r['b']:+.4f} (se {r['se']:.4f}) p={r['p']:.2g} "
              + (f"egger intercept={r['intercept']:+.4f} p={r['intercept_p']:.2g}" if "intercept" in r else ""))
    return out


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="eBMD:CAD,eBMD:MI,eBMD:STROKE,eBMD:HF,eBMD:OP")
    ap.add_argument("--reverse", default="CAD:eBMD,MI:eBMD,STROKE:eBMD,HF:eBMD,OP:eBMD")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()

    if a.self_test:
        self_test()
        print("\nwith directional pleiotropy (0.05):")
        self_test(pleiotropy=0.05, seed=11)
        return

    TAB.mkdir(parents=True, exist_ok=True)
    for label, spec in [("forward", a.pairs), ("reverse", a.reverse)]:
        rows, summ = [], []
        print(f"\n=== {label} MR")
        t0 = time.time()
        for pair in [p for p in spec.split(",") if p.strip()]:
            e, o = pair.split(":")
            r = run_pair(e.strip(), o.strip())
            if r is None:
                continue
            rows.append(r[0])
            summ.append(r[1])
        if rows:
            res = pd.concat(rows, ignore_index=True)
            res["or"] = np.where(res["outcome"].isin(BINARY), np.exp(res["b"]), np.nan)
            res["or_lo"] = np.where(res["outcome"].isin(BINARY), np.exp(res["b"] - 1.96 * res["se"]), np.nan)
            res["or_hi"] = np.where(res["outcome"].isin(BINARY), np.exp(res["b"] + 1.96 * res["se"]), np.nan)
            cols = ["exposure", "outcome", "method", "nsnp", "b", "se", "p", "or", "or_lo", "or_hi",
                    "q", "q_df", "q_p", "intercept", "intercept_se", "intercept_p"]
            res.reindex(columns=cols).to_csv(TAB / f"mr_results_{label}.csv", index=False, float_format="%.6g")
            pd.DataFrame(summ).to_csv(TAB / f"mr_sensitivity_{label}.csv", index=False, float_format="%.6g")
            print(f"  -> results/tables/02-mr/mr_results_{label}.csv  ({len(res)} rows, {time.time() - t0:.0f}s)")


if __name__ == "__main__":
    sys.exit(main())

