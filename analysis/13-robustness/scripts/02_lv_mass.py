"""Module 13B 鈥?does the eBMD effect reach the myocardium itself?

"Acts on the myocardium" is currently a conclusion by elimination: not coronary, not
circulating, not body size, not any annotated pathway. Left ventricular mass is the one
cardiac phenotype with a public GWAS big enough to test it directly.

  eBMD (592 instruments) -> left ventricular mass        GCST90244711  (n=41,855)
                          -> LV mass / body surface area GCST90435260  (n=36,083)

Both files are largely UK Biobank, so they overlap the eBMD exposure. Overlap inflates
false positives, which is why a positive result is reported with that caveat and a null is
the informative outcome.

Decision rule, fixed before the run:
  eBMD lowers LV mass with a CI excluding zero
        -> cardiac remodelling becomes measured evidence, not inference
  null with the CI excluding the equivalent of the heart-failure effect
        -> the axis is not concentric hypertrophy; it points at fibrosis, diastolic
           function or conduction, and the wording must say so

A secondary analysis, run only when the instrument table exists, asks whether LV mass
itself moves heart failure: LV mass instruments -> FinnGen heart failure.

Outputs (results/tables/13-robustness/):
  lv_mass.csv

Usage:
  python analysis/13-robustness/scripts/02_lv_mass.py
  python analysis/13-robustness/scripts/02_lv_mass.py --self-test
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from _common import estimate, have, implied_sd, mr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / "results" / "tables" / "13-robustness"
EXPOSURE = "eBMD"
OUTCOMES = {"LVM": "left ventricular mass (UKB CMR, grams)",
            "LVMI": "LV mass indexed to body surface area"}
REFERENCE_SD = 0.04          # a 0.04-SD shift in LV mass, the size the data suggest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    inst = mr.load_instruments(EXPOSURE)
    rows = []
    for key, label in OUTCOMES.items():
        r = estimate(EXPOSURE, key, inst=inst, label=f"eBMD -> {label}",
                     reference=REFERENCE_SD)
        if not r:
            continue
        sd = implied_sd(key)
        r["outcome_sd"] = sd
        r["b_sd"] = r["b"] / sd
        r["se_sd"] = r["se"] / sd
        r["lo_sd"] = r["lo"] / sd
        r["hi_sd"] = r["hi"] / sd
        from _power import mde, power
        r["power_ref_sd"] = power(r["se_sd"], REFERENCE_SD)
        r["mde80_sd"] = mde(r["se_sd"], 0.80)
        print(f"    standardised: {r['b_sd']:+.3f} SD per SD eBMD "
              f"({r['lo_sd']:+.3f}, {r['hi_sd']:+.3f}); implied trait SD {sd:.2f}; "
              f"power for {REFERENCE_SD} SD = {r['power_ref_sd']:.0%}")
        rows.append(r)
    if rows:
        pd.DataFrame(rows).to_csv(OUT / "lv_mass.csv", index=False, float_format="%.6g")
        print(f"\n-> {(OUT / 'lv_mass.csv').relative_to(ROOT).as_posix()}")
    return 0


def selftest():
    """A null cardiac outcome must stay null at a reported, adequate power.

    Scales match the real data: instrument effects on eBMD are ~0.04 with SE ~0.003, and
    a standardised cardiac trait has outcome SE ~0.02. That produces an IVW SE of ~0.02,
    which is what the module reports.
    """
    from _power import mde, power
    rng = np.random.default_rng(21)
    n = 592
    bx = rng.normal(0.04, 0.02, n)
    sx = np.full(n, 0.003)
    sy = np.full(n, 0.02)
    by = rng.normal(0, sy)
    r = mr.ivw(bx, by, sx, sy, random=True)
    p = power(r["se"], 0.05)
    print(f"selftest: null cardiac outcome b={r['b']:+.4f} se={r['se']:.4f} "
          f"p={r['p']:.3g}; power for 0.05 SD = {p:.0%}; MDE80 = {mde(r['se']):.3f}")
    ok = r["p"] > 0.01 and 0.3 < p < 1.0
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

