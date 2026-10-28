"""Follow-up on the one surviving causal signal: genetically predicted bone density
(eBMD) -> heart failure.

The module-02 analysis left eBMD -> HF as the only estimate that survived outlier
removal. This script asks whether that survives four further, largely independent
attacks:

  1. REPLICATION in two heart-failure GWAS that are not the discovery set:
       HF_SHAH     Shah 2020      (47,309 cases / 930,014 controls)
       HF_CHRONIC  Sakaue 2021    (14,262 cases / 471,898 controls)
  2. STEIGER-FILTERED re-estimation: drop variants that explain more variance in HF
     than in eBMD, then re-fit IVW and MR-Egger.
  3. MULTIVARIABLE MR: does the eBMD effect survive adjustment for body mass index,
     systolic blood pressure, type-2 diabetes and coronary artery disease? A direct
     effect that survives is a stronger claim; an effect that vanishes identifies the
     route (adiposity / vascular / glycaemic).
  4. SHARED-LOCUS check: remove eBMD instruments that sit within 500 kb of an HF
     instrument, i.e. loci already known to influence HF, and re-estimate.

Outputs (results/tables/05-hf-followup/):
  hf_replication.csv
  hf_steiger.csv
  hf_mvmr.csv
  hf_shared_loci.csv
  hf_outlier_loci.csv

Usage:
  python analysis/05-hf-followup/scripts/01_hf_followup.py
  python analysis/05-hf-followup/scripts/01_hf_followup.py --self-test
"""
import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
MR_DIR = HERE.parents[1] / "02-mr" / "scripts"
sys.path.insert(0, str(MR_DIR))

_spec = importlib.util.spec_from_file_location("mr", MR_DIR / "03_mr.py")
mr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mr)

ROOT = HERE.parents[2]
OUT = ROOT / "results" / "tables" / "05-hf-followup"
EXPOSURE = "eBMD"
OUTCOME = "HF"
REPLICATES = ["HF_SHAH", "HF_CHRONIC"]
MEDIATORS = ["BMI", "SBP", "T2D", "CAD"]


# ---------------------------------------------------------------- helpers
def align(ea, oa, beta, ref_ea, ref_oa):
    """Return beta oriented to (ref_ea, ref_oa), or None if the alleles are incompatible."""
    if (ea, oa) == (ref_ea, ref_oa):
        return beta
    if (ea, oa) == (ref_oa, ref_ea):
        return -beta
    return None


def stack(exposure_key, outcome_key, extra_keys=()):
    """Build one table where every effect is oriented to the exposure effect allele."""
    inst = mr.load_instruments(exposure_key)
    if inst is None or inst.empty:
        return None
    snps = inst["SNP"].tolist()
    out = mr.outcome_at_snps(outcome_key, snps)
    if out is None or out.empty:
        return None
    base = mr.harmonise(inst, out)                      # eBMD vs outcome, aligned
    if base.empty:
        return None
    keep = base.set_index("SNP")[["effect_allele", "other_allele", "bx", "sx", "by", "sy",
                                  "n_x", "n_y", "px"]]
    for k in extra_keys:
        tab = mr.outcome_at_snps(k, snps)
        if tab is None or tab.empty:
            keep[f"b_{k}"] = np.nan
            keep[f"s_{k}"] = np.nan
            continue
        rows = []
        for snp, r in keep.iterrows():
            t = tab[tab["SNP"] == snp]
            if t.empty:
                rows.append((snp, np.nan, np.nan))
                continue
            t = t.iloc[0]
            b = align(t["effect_allele"], t["other_allele"], float(t["beta"]),
                      r["effect_allele"], r["other_allele"])
            rows.append((snp, b, float(t["se"]) if b is not None else np.nan))
        d = pd.DataFrame(rows, columns=["SNP", f"b_{k}", f"s_{k}"]).set_index("SNP")
        keep = keep.join(d)
    return keep.reset_index()


def ivw_from(bx, by, sx, sy):
    return mr.ivw(bx, by, sx, sy, random=True)


def steiger_filter(tab):
    """Per-SNP variance explained in exposure vs outcome; keep exposure-dominant SNPs."""
    def r2(b, se, n):
        t = b / se
        r = t / np.sqrt(t ** 2 + np.maximum(n - 2, 1))
        return r ** 2
    rx = r2(tab["bx"].to_numpy(float), tab["sx"].to_numpy(float), tab["n_x"].to_numpy(float))
    ry = r2(tab["by"].to_numpy(float), tab["sy"].to_numpy(float), tab["n_y"].to_numpy(float))
    return rx > ry, rx, ry


# ---------------------------------------------------------------- MVMR
def mvmr_ivw(X, y, sy):
    """Multivariable IVW with sandwich standard errors.

    X: n x k matrix of exposure effects on the instruments; y: outcome effects;
    sy: outcome standard errors. Returns coefficients and SEs.
    """
    w = 1.0 / sy ** 2
    W = np.diag(w)
    XtWX = X.T @ W @ X
    if np.linalg.matrix_rank(XtWX) < X.shape[1]:
        return None
    beta = np.linalg.solve(XtWX, X.T @ W @ y)
    resid = y - X @ beta
    meat = X.T @ W @ np.diag(resid ** 2) @ W @ X
    cov = np.linalg.solve(XtWX, meat @ np.linalg.solve(XtWX, XtWX).T) if False else \
        np.linalg.inv(XtWX) @ meat @ np.linalg.inv(XtWX)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    return beta, se, resid


def conditional_f(bx1, sx1, bx2, n):
    """Approximate conditional F for exposure 1 given exposure 2 (Sanderson-style):
    the part of bx1 not explained by bx2, scaled by its own precision."""
    v = np.column_stack([np.ones_like(bx2), bx2])
    wts = 1.0 / sx1 ** 2
    Wv = v * wts[:, None]
    coef = np.linalg.lstsq(Wv.T @ v, Wv.T @ bx1, rcond=None)[0]
    resid = bx1 - v @ coef
    q = np.sum(wts * resid ** 2)
    return float(q / max(n - 2, 1))


# ---------------------------------------------------------------- analyses
def replication():
    rows = []
    for key in [OUTCOME] + REPLICATES:
        tab = stack(EXPOSURE, key)
        if tab is None:
            print(f"  {key}: no data")
            continue
        bx, by = tab["bx"].to_numpy(float), tab["by"].to_numpy(float)
        sx, sy = tab["sx"].to_numpy(float), tab["sy"].to_numpy(float)
        r = ivw_from(bx, by, sx, sy)
        e = mr.egger(bx, by, sx, sy)
        rows.append(dict(outcome=key, nsnp=len(tab), b=r["b"], se=r["se"], p=r["p"],
                         lo=r["b"] - 1.96 * r["se"], hi=r["b"] + 1.96 * r["se"],
                         or_=np.exp(r["b"]), or_lo=np.exp(r["b"] - 1.96 * r["se"]),
                         or_hi=np.exp(r["b"] + 1.96 * r["se"]),
                         q_p=r["q_p"], egger_b=e["b"], egger_p=e["p"],
                         egger_intercept_p=e["intercept_p"]))
        print(f"  {key:11s} n={len(tab):4d}  b={r['b']:+.4f} (p={r['p']:.2g})  "
              f"OR={np.exp(r['b']):.3f}  Egger intercept p={e['intercept_p']:.2g}")
    return pd.DataFrame(rows)


def steiger(tab):
    keep, rx, ry = steiger_filter(tab)
    rows = []
    for label, mask in [("all variants", np.ones(len(tab), bool)), ("Steiger-passing", keep)]:
        bx, by = tab["bx"].to_numpy(float)[mask], tab["by"].to_numpy(float)[mask]
        sx, sy = tab["sx"].to_numpy(float)[mask], tab["sy"].to_numpy(float)[mask]
        if mask.sum() < 5:
            continue
        r = ivw_from(bx, by, sx, sy)
        e = mr.egger(bx, by, sx, sy)
        rows.append(dict(subset=label, nsnp=int(mask.sum()), b=r["b"], se=r["se"], p=r["p"],
                         q_p=r["q_p"], egger_p=e["p"], egger_intercept_p=e["intercept_p"]))
        print(f"  {label:16s} n={int(mask.sum()):4d}  b={r['b']:+.4f} (p={r['p']:.2g})")
    return pd.DataFrame(rows), keep


def mvmr(tab):
    rows = []
    base = ivw_from(tab["bx"].to_numpy(float), tab["by"].to_numpy(float),
                    tab["sx"].to_numpy(float), tab["sy"].to_numpy(float))
    for med in MEDIATORS:
        col = f"b_{med}"
        if col not in tab:
            print(f"  {med}: mediator data missing")
            continue
        sub = tab.dropna(subset=[col, f"s_{med}"])
        if len(sub) < 20:
            print(f"  {med}: only {len(sub)} shared instruments, skipped")
            continue
        X = np.column_stack([sub["bx"].to_numpy(float), sub[col].to_numpy(float)])
        y = sub["by"].to_numpy(float)
        sy = sub["sy"].to_numpy(float)
        res = mvmr_ivw(X, y, sy)
        if res is None:
            print(f"  {med}: singular design matrix")
            continue
        beta, se, resid = res
        p = 2 * stats.norm.sf(np.abs(beta / se))
        cf = conditional_f(sub["bx"].to_numpy(float), sub["sx"].to_numpy(float),
                           sub[col].to_numpy(float), len(sub))
        corr = float(np.corrcoef(sub["bx"].to_numpy(float), sub[col].to_numpy(float))[0, 1])
        rows.append(dict(mediator=med, nsnp=len(sub),
                         direct_b=beta[0], direct_se=se[0], direct_p=p[0],
                         direct_or=np.exp(beta[0]),
                         direct_or_lo=np.exp(beta[0] - 1.96 * se[0]),
                         direct_or_hi=np.exp(beta[0] + 1.96 * se[0]),
                         mediator_b=beta[1], mediator_se=se[1], mediator_p=p[1],
                         univariable_b=base["b"], univariable_p=base["p"],
                         attenuation_pct=100 * (1 - beta[0] / base["b"]),
                         conditional_F_eBMD=cf,
                         corr_effects=corr,
                         resid_scale=float(np.sqrt(np.mean(resid ** 2)))))
        print(f"  {med:4s} n={len(sub):4d}  eBMD direct b={beta[0]:+.4f} (p={p[0]:.2g}) "
              f"vs univariable {base['b']:+.4f}; {med} b={beta[1]:+.4f} (p={p[1]:.2g})")
    return pd.DataFrame(rows)


def shared_locus_test(tab):
    hf = mr.load_instruments(OUTCOME)
    if hf is None or hf.empty:
        return pd.DataFrame(), None
    inst = mr.load_instruments(EXPOSURE)
    pos = inst.set_index("SNP")[["chr", "pos"]]
    hfpos = hf.set_index("SNP")[["chr", "pos"]]
    dist = {}
    for snp in tab["SNP"]:
        if snp not in pos.index:
            dist[snp] = np.inf
            continue
        c, p = float(pos.loc[snp, "chr"]), float(pos.loc[snp, "pos"])
        same = hfpos[hfpos["chr"].astype(float) == c]
        dist[snp] = float(np.min(np.abs(same["pos"].astype(float) - p))) if len(same) else np.inf
    tab = tab.copy()
    tab["dist_to_hf_instrument"] = tab["SNP"].map(dist)
    print(f"  eBMD instruments within 500 kb of an HF instrument: "
          f"{int((tab['dist_to_hf_instrument'] <= 500_000).sum())} of {len(tab)}")
    rows = []
    for label, mask in [("all", np.ones(len(tab), bool)),
                        ("excluding HF-near loci", tab["dist_to_hf_instrument"] > 500_000)]:
        bx, by = tab["bx"].to_numpy(float)[mask], tab["by"].to_numpy(float)[mask]
        sx, sy = tab["sx"].to_numpy(float)[mask], tab["sy"].to_numpy(float)[mask]
        if mask.sum() < 5:
            continue
        r = ivw_from(bx, by, sx, sy)
        rows.append(dict(subset=label, nsnp=int(mask.sum()), b=r["b"], se=r["se"], p=r["p"],
                         q_p=r["q_p"]))
        print(f"  {label:24s} n={int(mask.sum()):4d}  b={r['b']:+.4f} (p={r['p']:.2g})")
    return pd.DataFrame(rows), tab


def outlier_loci(tab, max_show=8):
    keep, removed = presso_mask(tab)
    if not removed:
        return pd.DataFrame()
    sub = tab.iloc[removed].copy()
    sub["z"] = (sub["by"] / sub["bx"])  # descriptive
    return sub[["SNP", "by", "bx"]].head(max_show)


def presso_mask(tab):
    bx, by = tab["bx"].to_numpy(float), tab["by"].to_numpy(float)
    sx, sy = tab["sx"].to_numpy(float), tab["sy"].to_numpy(float)
    import importlib.util as iu
    spec = iu.spec_from_file_location("pl", MR_DIR / "04_pleiotropy.py")
    pl = iu.module_from_spec(spec)
    spec.loader.exec_module(pl)
    return pl.presso_outliers(bx, by, sx, sy)


# ---------------------------------------------------------------- self test
def selftest():
    """MVMR must recover a known direct effect and the conditional F must be finite."""
    rng = np.random.default_rng(11)
    n = 300
    b1 = rng.normal(0, 0.05, n)
    b2 = 0.35 * b1 + rng.normal(0, 0.05, n)          # correlated exposures
    s1 = np.full(n, 0.004)
    s2 = np.full(n, 0.006)
    sy = np.full(n, 0.02)
    y = 0.25 * b1 + 0.40 * b2 + rng.normal(0, sy)    # true direct effects
    beta, se, resid = mvmr_ivw(np.column_stack([b1, b2]), y, sy)
    cf = conditional_f(b1, s1, b2, n)
    print(f"selftest: MVMR direct eBMD {beta[0]:+.3f} (se {se[0]:.3f}) vs truth 0.25")
    print(f"selftest: MVMR mediator {beta[1]:+.3f} (se {se[1]:.3f}) vs truth 0.40")
    print(f"selftest: conditional F = {cf:.1f}")
    ok = abs(beta[0] - 0.25) < 0.08 and abs(beta[1] - 0.40) < 0.10 and np.isfinite(cf) and cf > 0
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    print("== 1. replication in independent heart-failure GWAS")
    rep = replication()
    rep.to_csv(OUT / "hf_replication.csv", index=False, float_format="%.6g")

    # the MVMR/Steiger/shared-locus analyses all need the mediator effects on the same
    # instrument set, so they are stacked once here
    tab = stack(EXPOSURE, OUTCOME, MEDIATORS)
    print(f"\nstacked {len(tab)} instruments with {len(tab.columns)} columns: "
          f"{', '.join(c for c in tab.columns if c.startswith('b_'))}")

    print("\n== 2. Steiger-filtered re-estimation")
    st, keep = steiger(tab)
    st.to_csv(OUT / "hf_steiger.csv", index=False, float_format="%.6g")

    print("\n== 3. multivariable MR (adjusting for candidate mediators)")
    mv = mvmr(tab)
    if len(mv):
        mv.to_csv(OUT / "hf_mvmr.csv", index=False, float_format="%.6g")

    print("\n== 4. shared-locus check")
    sl, tab2 = shared_locus_test(tab)
    if len(sl):
        sl.to_csv(OUT / "hf_shared_loci.csv", index=False, float_format="%.6g")

    print("\n== outliers removed by the MR-PRESSO-style step")
    ol = outlier_loci(tab)
    if len(ol):
        ol.to_csv(OUT / "hf_outlier_loci.csv", index=False, float_format="%.6g")
        print(ol.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())

