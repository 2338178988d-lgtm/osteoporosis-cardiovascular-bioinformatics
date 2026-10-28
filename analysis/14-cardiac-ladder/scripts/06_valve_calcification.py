"""Module 14, part 4 鈥?the valve question: calcific or degenerative?

FinnGen publishes only three valve endpoints (I9_VHD 28,602; I9_CAVS_OPERATED 12,418;
I9_RHEUVALV 877) and no mitral / aortic-regurgitation codes, so the protective estimate on
"I9_VHD" cannot be split by lesion type from Finnish data alone. That matters:

  gain-of-mineralisation reading  -> calcific disease drives it; calcific aortic stenosis
                                     should be protective, and calcific aortic valve
                                     calcification, measured directly, more so
  matrix / degeneration reading   -> the calcific subset is only part of valvular disease;
                                     direct valve calcification may not move, or may move less

Module 14's indirect evidence was ambiguous: I9_CAVS_OPERATED -0.024 (p=0.47, 47 % power),
too underpowered to decide. GCST90310293 measures aortic valve calcification directly in
14,819 cases / 941,863 controls (PMID 38494474), which is a larger and cleaner case
definition than "operated CAVS".

DECISION RULE, fixed before the run:
  AVC protective at p < 0.05
     -> the valvular effect is calcific; the gain-of-mineralisation reading survives and the
        module-13 "detect ectopic calcification" prediction is supported
  AVC null with a CI that excludes -0.074 (the I9_VHD estimate)
     -> the valvular effect is NOT calcific; I9_VHD must be read as degenerative/mixed
        disease, and the calcification reading is confined to whatever remains
  AVC null with a CI containing -0.074
     -> underpowered, report the power and leave the question open

Exposure: the same 600 eBMD instruments. The AVC file carries an `n` column whose median is
about 407,000 (a per-variant analysis sample size), so instrument R^2 is computed from the
file's own `n` where available.

Outputs (results/tables/14-cardiac-ladder/):
  avc.csv     the eBMD -> AVC estimate with power, next to the FinnGen valve rungs

Usage:
  python analysis/14-cardiac-ladder/scripts/06_valve_calcification.py
  python analysis/14-cardiac-ladder/scripts/06_valve_calcification.py --self-test
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "analysis" / "13-robustness" / "scripts"))

from _common import estimate, have, mr  # noqa: E402

OUT = ROOT / "results" / "tables" / "14-cardiac-ladder"
EXPOSURE = "eBMD"
COMPARATORS = {"FG_VHD": "valvular heart disease, ICD aggregate",
               "FG_CAVS": "calcific aortic stenosis, operated"}


def decide(s, reference=-0.074):
    """Classify one estimate against the pre-registered rule.

    "CI excludes the I9_VHD estimate" means the interval lies wholly above or wholly below
    that value 鈥?for a null estimate the informative case is wholly above it.
    """
    if s["p"] < 0.05 and s["b"] < 0:
        return "protective"
    if s["p"] >= 0.05 and (s["lo"] > reference or s["hi"] < reference):
        return "not calcific"
    return "open"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    inst = mr.load_instruments(EXPOSURE)
    if inst is None or inst.empty:
        print(f"{EXPOSURE}: no instruments")
        return 2
    if not have("AVC"):
        print("AVC: file absent 鈥?run download_plan.py --only AVC")
        return 2

    rows = []
    r = estimate(EXPOSURE, "AVC", inst=inst, label="heel eBMD -> aortic valve calcification")
    if r:
        r["role"] = "direct calcification endpoint"
        r["cases"] = 14819
        rows.append(r)

    main_row = r
    if main_row:
        print("\n--- decision rule")
        verdict = decide(main_row)
        if verdict == "protective":
            print("  AVC protective: the valvular effect is calcific. The gain-of-"
                  "mineralisation reading survives, and module 13's prediction that eBMD "
                  "should touch ectopic calcification is supported.")
        elif verdict == "not calcific":
            print("  AVC null with the CI excluding -0.074: the valvular effect is NOT "
                  "calcific. I9_VHD must be read as degenerative/mixed disease.")
        else:
            print(f"  AVC null and underpowered ({main_row['power_ref']:.0%} power), CI "
                  "spanning -0.074. The question stays open; report it as insufficient "
                  "rather than negative.")

    # the two Finnish valve rungs, for side-by-side reading
    ladder = OUT / "ladder.csv"
    if ladder.exists():
        lad = pd.read_csv(ladder).set_index("outcome")
        for key, label in COMPARATORS.items():
            if key in lad.index:
                s = lad.loc[key]
                rows.append(dict(pair=f"{EXPOSURE} -> {key}", exposure=EXPOSURE,
                                 outcome=key, nsnp=int(s["nsnp"]), b=float(s["b"]),
                                 se=float(s["se"]), p=float(s["p"]), lo=float(s["lo"]),
                                 hi=float(s["hi"]), q_p=float(s["q_p"]),
                                 egger_intercept_p=float(s["egger_intercept_p"]),
                                 power_ref=float(s["power_ref"]),
                                 mde_80=float(s["mde_80"]),
                                 role=label, cases=int(s["cases"])))

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "avc.csv", index=False, float_format="%.6g")
    print(f"\n-> {(OUT / 'avc.csv').relative_to(ROOT).as_posix()}")
    print("\n  valve endpoints side by side")
    for _, s in df.iterrows():
        print(f"    {s['role'][:44]:44s} b={s['b']:+.4f} "
              f"({s['lo']:+.4f}, {s['hi']:+.4f}) p={s['p']:.3g} "
              f"power={s['power_ref']:.0%}")
    return 0


def selftest():
    """The decision rule must classify three planted cases correctly."""
    cases = [("protective", dict(b=-0.08, se=0.03, p=0.008, lo=-0.14, hi=-0.02,
                                power_ref=0.6), "protective"),
             ("null, CI excludes -0.074", dict(b=0.01, se=0.025, p=0.69, lo=-0.04,
                                               hi=0.06, power_ref=0.9), "not calcific"),
             ("null, wide CI", dict(b=-0.02, se=0.06, p=0.74, lo=-0.14, hi=0.10,
                                    power_ref=0.3), "open")]
    ok = True
    for name, s, expect in cases:
        got = decide(s)
        good = got == expect
        ok &= good
        print(f"selftest {name:26s} -> {got:14s} {'ok' if good else 'WRONG'}")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

