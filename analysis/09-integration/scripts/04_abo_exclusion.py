"""Module 09 step 4 鈥?is the bone-density effect on heart failure carried by ABO?

The colocalisation scan found exactly one locus out of 600 where eBMD and heart failure
share a causal variant, and it is ABO. ABO is the most pleiotropic locus in the human
genome: it sets blood-group glycosylation and shifts von Willebrand factor, factor VIII,
soluble E-selectin, ICAM-1, LDL and dozens of other traits. A shared causal variant there
is therefore at least as likely to be horizontal pleiotropy as a bone-to-heart mechanism.

That makes one test decisive rather than decorative: **remove every instrument within a
megabase of ABO and re-estimate**. If the heart-failure effect collapses, the project's
one surviving causal claim is ABO pleiotropy and must be reported as such. If it survives,
ABO is one contributor among many.

Also tested: removal of the six loci where both traits associate but with *different*
causal variants (PP.H3 >= 0.8), which are by construction not shared mechanisms.

Outputs (results/tables/09-integration/):
  abo_exclusion.csv     the effect before and after each exclusion, both cohorts

Usage:
  python analysis/09-integration/scripts/04_abo_exclusion.py
  python analysis/09-integration/scripts/04_abo_exclusion.py --self-test
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
OUT = ROOT / "results" / "tables" / "09-integration"
EXPOSURE = "eBMD"
OUTCOMES = [("HF", "HERMES (overlapping)"), ("FG_HF", "FinnGen R12 (independent)")]
ABO = (9, 133_274_306)
ABO_WINDOW = 1_000_000
H3_LOCI = [("rs11657636", 17, 2150626), ("rs67631072", 1, 37996149), ("rs4129645", 5, 42526795),
           ("rs7197197", 16, 72877507), ("rs115108083", 1, 50656577), ("rs80195615", 1, 61180831)]


def estimate(d):
    bx, by = d["bx"].to_numpy(float), d["by"].to_numpy(float)
    sx, sy = d["sx"].to_numpy(float), d["sy"].to_numpy(float)
    r = mr.ivw(bx, by, sx, sy, random=True)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    inst = mr.load_instruments(EXPOSURE)
    inst["chr"] = inst["chr"].astype(float)
    inst["pos"] = inst["pos"].astype(float)

    near_abo = (inst["chr"] == ABO[0]) & (inst["pos"].between(ABO[1] - ABO_WINDOW,
                                                              ABO[1] + ABO_WINDOW))
    near_h3 = pd.Series(False, index=inst.index)
    for _, c, p in H3_LOCI:
        near_h3 |= (inst["chr"] == c) & (inst["pos"].between(p - 500_000, p + 500_000))
    print(f"instruments within 1 Mb of ABO: {int(near_abo.sum())}")
    print(f"instruments within 500 kb of an H3 locus: {int(near_h3.sum())}")

    rows = []
    for key, label in OUTCOMES:
        out = mr.outcome_at_snps(key, inst["SNP"].tolist())
        if out is None or out.empty:
            print(f"  {key}: no outcome data")
            continue
        d = mr.harmonise(inst, out)
        d["chr"] = d["SNP"].map(dict(zip(inst["SNP"], inst["chr"])))
        d["pos"] = d["SNP"].map(dict(zip(inst["SNP"], inst["pos"])))
        in_abo = (d["chr"] == ABO[0]) & (d["pos"].between(ABO[1] - ABO_WINDOW, ABO[1] + ABO_WINDOW))
        in_h3 = pd.Series(False, index=d.index)
        for _, c, p in H3_LOCI:
            in_h3 |= (d["chr"] == c) & (d["pos"].between(p - 500_000, p + 500_000))

        for subset, mask in [("all instruments", pd.Series(True, index=d.index)),
                             ("ABO region removed", ~in_abo),
                             ("H3 loci removed", ~in_h3),
                             ("ABO + H3 removed", ~in_abo & ~in_h3)]:
            sub = d[mask]
            if len(sub) < 20:
                continue
            r = estimate(sub)
            rows.append(dict(cohort=label, outcome=key, subset=subset, nsnp=len(sub),
                             n_removed=int((~mask).sum()),
                             b=r["b"], se=r["se"], p=r["p"],
                             lo=r["b"] - 1.96 * r["se"], hi=r["b"] + 1.96 * r["se"],
                             or_=float(np.exp(r["b"])), q_p=r["q_p"]))
            print(f"  {label:28s} {subset:22s} n={len(sub):4d} "
                  f"b={r['b']:+.4f} ({r['b'] - 1.96 * r['se']:+.4f}, "
                  f"{r['b'] + 1.96 * r['se']:+.4f})  p={r['p']:.3g}")

    d = pd.DataFrame(rows)
    if len(d):
        # change caused by the exclusion, relative to the all-instrument estimate
        d["change_vs_all"] = np.nan
        for (coh, out), g in d.groupby(["cohort", "outcome"]):
            base = g.loc[g.subset == "all instruments", "b"]
            if len(base):
                d.loc[g.index, "change_vs_all"] = g["b"] - float(base.iloc[0])
        d.to_csv(OUT / "abo_exclusion.csv", index=False, float_format="%.6g")
        print(f"\n-> {OUT / 'abo_exclusion.csv'}")
    return 0


def selftest():
    """Removing a set of instruments with a large common effect must move the estimate;
    removing null instruments must not."""
    rng = np.random.default_rng(5)
    n = 300
    bx = rng.normal(0, 0.05, n)
    sx = np.full(n, 0.004)
    sy = np.full(n, 0.02)
    by = -0.03 * bx + rng.normal(0, sy)          # true effect -0.03
    d = pd.DataFrame(dict(bx=bx, sx=sx, by=by, sy=sy))
    base = estimate(d)["b"]
    # replace 20 instruments with a strong spurious effect
    d2 = d.copy()
    d2.loc[:19, "by"] = 0.35 * d2.loc[:19, "bx"]
    with_spurious = estimate(d2)["b"]
    trimmed = estimate(d2.iloc[20:])["b"]
    print(f"selftest: clean {base:+.4f}, with spurious {with_spurious:+.4f}, "
          f"after removing them {trimmed:+.4f}")
    ok = abs(trimmed - base) < 0.01 and abs(with_spurious) < abs(base) / 2
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

