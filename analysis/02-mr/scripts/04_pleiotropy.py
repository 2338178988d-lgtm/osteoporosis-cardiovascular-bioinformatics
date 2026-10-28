"""Pleiotropy sensitivity for the 02-mr estimates.

The LD-clumped IVW estimates for eBMD -> CAD / HF come with a significant MR-Egger
intercept, which means directional pleiotropy and makes the IVW estimate unreliable.
This script tests whether the effect survives after removing the variants that carry
that pleiotropy, and how much of it depends on the strongest loci.

Methods
  MR-PRESSO-style outlier removal: regress by on bx, compute per-SNP residuals, and
  iteratively drop variants whose residual exceeds the Bonferroni-corrected threshold
  (the same criterion MR-PRESSO uses), then re-estimate with IVW.
  Cook's-distance style leave-one-locus-out: re-estimate after removing each of the 10
  strongest loci in turn, to show the estimate is not driven by a single region.
  Egger intercept on the trimmed set, to check the pleiotropy is actually gone.

Outputs (results/tables/02-mr/):
  mr_pleiotropy_sensitivity.csv   per pair: full / trimmed IVW, outliers removed, Egger
  mr_loco_influence.csv           leave-one-locus-out estimates for the flagged pairs

Usage:
    python analysis/02-mr/scripts/04_pleiotropy.py
    python analysis/02-mr/scripts/04_pleiotropy.py --self-test
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import _sumstats
import importlib.util

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("mr", HERE / "03_mr.py")
mr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mr)

TAB = mr.TAB


def presso_outliers(bx, by, sx, sy, max_iter=10):
    """Iterative residual outlier removal; returns the mask of variants kept."""
    keep = np.ones(len(bx), bool)
    removed = []
    for _ in range(max_iter):
        b, s = bx[keep], by[keep]
        se = sy[keep] / np.abs(b)
        w = 1 / se ** 2
        est = np.sum(w * (s / b)) / np.sum(w)
        resid = (s / b) - est
        # the outlier test uses the ratio SE, as MR-PRESSO does on the per-SNP estimates
        z = resid / se
        # Bonferroni over the variants currently in the set
        thr = stats.norm.isf(0.025 / max(len(z), 1))
        bad = np.abs(z) > thr
        if not bad.any():
            break
        idx = np.where(keep)[0][bad]
        keep[idx] = False
        removed.extend(idx.tolist())
    return keep, removed


def analyse(exp_key, out_key):
    inst = mr.load_instruments(exp_key)
    if inst is None or inst.empty:
        return None
    out = mr.outcome_at_snps(out_key, inst["SNP"].tolist())
    if out is None or out.empty:
        return None
    d = mr.harmonise(inst, out)
    if len(d) < 5:
        return None
    bx, by = d["bx"].to_numpy(), d["by"].to_numpy()
    sx, sy = d["sx"].to_numpy(), d["sy"].to_numpy()

    full = mr.ivw(bx, by, sx, sy)
    egg_full = mr.egger(bx, by, sx, sy)
    keep, removed = presso_outliers(bx, by, sx, sy)
    trim = mr.ivw(bx[keep], by[keep], sx[keep], sy[keep]) if keep.sum() >= 5 else None
    egg_trim = mr.egger(bx[keep], by[keep], sx[keep], sy[keep]) if keep.sum() >= 5 else None

    row = dict(exposure=exp_key, outcome=out_key, nsnp=len(d),
               n_outliers=len(removed),
               outlier_snps=";".join(d["SNP"].to_numpy()[removed][:20]),
               ivw_b=full["b"], ivw_se=full["se"], ivw_p=full["p"],
               ivw_q_p=full["q_p"],
               egger_b=egg_full["b"], egger_p=egg_full["p"],
               egger_intercept=egg_full["intercept"], egger_intercept_p=egg_full["intercept_p"])
    if trim is not None:
        row.update(trim_nsnp=int(keep.sum()), trim_ivw_b=trim["b"], trim_ivw_se=trim["se"],
                   trim_ivw_p=trim["p"], trim_q_p=trim["q_p"],
                   trim_egger_intercept=egg_trim["intercept"],
                   trim_egger_intercept_p=egg_trim["intercept_p"])
        if mr.BINARY and out_key in mr.BINARY:
            row["trim_or"] = float(np.exp(trim["b"]))
            row["trim_or_lo"] = float(np.exp(trim["b"] - 1.96 * trim["se"]))
            row["trim_or_hi"] = float(np.exp(trim["b"] + 1.96 * trim["se"]))
    if out_key in mr.BINARY:
        row["ivw_or"] = float(np.exp(full["b"]))
        row["ivw_or_lo"] = float(np.exp(full["b"] - 1.96 * full["se"]))
        row["ivw_or_hi"] = float(np.exp(full["b"] + 1.96 * full["se"]))
    return row, d, keep, removed, inst


def loco(d, keep, inst, top_n=10):
    """Leave-one-locus-out: drop each of the top_n strongest variants and re-estimate."""
    bx, by = d["bx"].to_numpy(), d["by"].to_numpy()
    sx, sy = d["sx"].to_numpy(), d["sy"].to_numpy()
    order = np.argsort(d["px"].to_numpy())
    loc = inst.set_index("SNP")
    base = mr.ivw(bx[keep], by[keep], sx[keep], sy[keep])
    rows = []
    for i in order[:top_n]:
        m = keep.copy()
        m[i] = False
        if m.sum() < 5:
            continue
        r = mr.ivw(bx[m], by[m], sx[m], sy[m])
        snp = d["SNP"].to_numpy()[i]
        rows.append(dict(removed_snp=snp,
                         removed_chr=int(loc.loc[snp, "chr"]) if snp in loc.index else -1,
                         removed_pos=int(loc.loc[snp, "pos"]) if snp in loc.index else -1,
                         base_b=base["b"], loco_b=r["b"],
                         shift=r["b"] - base["b"], loco_p=r["p"]))
    return pd.DataFrame(rows)


PAIRS = [("eBMD", "CAD"), ("eBMD", "HF"), ("eBMD", "MI"), ("eBMD", "STROKE"),
         ("eBMD", "OP"), ("OP", "CAD"), ("OP", "MI"), ("OP", "STROKE"), ("OP", "HF"),
         ("CAD", "eBMD"), ("MI", "eBMD"), ("STROKE", "eBMD"), ("HF", "eBMD"), ("OP", "eBMD"),
         ("CAD", "OP"), ("MI", "OP"), ("STROKE", "OP"), ("HF", "OP")]


def selftest():
    """A causal effect plus a few gross outlier variants: PRESSO must find them."""
    rng = np.random.default_rng(5)
    n = 200
    bx = rng.uniform(0.02, 0.10, n) * rng.choice([-1, 1], n)
    sx = np.full(n, 0.003)
    sy = np.full(n, 0.02)
    by = 0.30 * bx + rng.normal(0, sy)
    bad = rng.choice(n, 6, replace=False)
    by[bad] += rng.choice([-1, 1], 6) * 0.09          # gross pleiotropy on 6 variants
    keep, removed = presso_outliers(bx, by, sx, sy)
    found = len(set(removed) & set(bad))
    full = mr.ivw(bx, by, sx, sy)["b"]
    trim = mr.ivw(bx[keep], by[keep], sx[keep], sy[keep])["b"]
    print(f"selftest: {found}/6 planted outliers found, {removed.count(0) if 0 in removed else len(removed)} removed total")
    print(f"selftest: full IVW b={full:.3f} -> trimmed b={trim:.3f} (truth 0.30)")
    ok = found >= 4 and abs(trim - 0.30) < 0.02
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    rows, loco_frames = [], []
    for e, o in PAIRS:
        r = analyse(e, o)
        if r is None:
            print(f"{e} -> {o}: skipped (missing data)")
            continue
        row, d, keep, removed, inst = r
        rows.append(row)
        print(f"{e} -> {o}: n={row['nsnp']}, outliers={row['n_outliers']}, "
              f"IVW {row['ivw_b']:+.3f} (p={row['ivw_p']:.2g}) -> "
              f"trimmed {row.get('trim_ivw_b', float('nan')):+.3f} "
              f"(p={row.get('trim_ivw_p', float('nan')):.2g}); "
              f"Egger intercept p {row['egger_intercept_p']:.2g} -> "
              f"{row.get('trim_egger_intercept_p', float('nan')):.2g}")
        if row["egger_intercept_p"] < 0.05 or row["ivw_p"] < 0.05:
            l = loco(d, keep, inst)
            if len(l):
                l.insert(0, "outcome", o)
                l.insert(0, "exposure", e)
                loco_frames.append(l)

    s = pd.DataFrame(rows)
    s.to_csv(TAB / "mr_pleiotropy_sensitivity.csv", index=False, float_format="%.6g")
    print(f"\n-> results/tables/02-mr/mr_pleiotropy_sensitivity.csv ({len(s)} pairs)")
    if loco_frames:
        lo = pd.concat(loco_frames, ignore_index=True)
        lo.to_csv(TAB / "mr_loco_influence.csv", index=False, float_format="%.6g")
        print(f"-> results/tables/02-mr/mr_loco_influence.csv ({len(lo)} rows)")
        print("\nlargest leave-one-locus-out shifts:")
        print(lo.reindex(lo["shift"].abs().sort_values(ascending=False).index)
              .head(10)[["exposure", "outcome", "removed_snp", "removed_chr",
                         "base_b", "loco_b", "shift"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())

