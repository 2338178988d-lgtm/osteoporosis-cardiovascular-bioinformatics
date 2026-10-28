"""Module 14, one number instead of three: does the ladder behave like a dilution series?

If heart failure is a mixture of an ischaemic part (where the effect is absent) and a
non-ischaemic part (where it is present), then the all-HF estimate and the HF+CHD estimate
together imply what the non-ischaemic part must be. That implied value can be compared with
the separately measured non-ischaemic cardiomyopathy estimate. Agreement means the three
numbers are one number.

Output: results/tables/14-cardiac-ladder/dilution.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
LAD = ROOT / "results" / "tables" / "14-cardiac-ladder" / "ladder.csv"
SUB = ROOT / "results" / "tables" / "05-hf-followup" / "fg_subtype.csv"
OUT = ROOT / "results" / "tables" / "14-cardiac-ladder" / "dilution.csv"

# risk sets are not full partitions, so cases are used only as the mixing weight
CASES_ALL_HF = 37_653
CASES_HF_CHD = 17_532


def main():
    lad = pd.read_csv(LAD).set_index("outcome")
    sub = pd.read_csv(SUB).set_index("outcome")
    hf, chd = sub.loc["FG_HF"], sub.loc["FG_HF_CHD"]
    w = CASES_HF_CHD / CASES_ALL_HF

    # all-HF = w * HF+CHD + (1 - w) * HF-without-CHD
    b_non = (hf["b"] - w * chd["b"]) / (1 - w)
    se_non = float(np.sqrt(hf["se"] ** 2 + (w * chd["se"]) ** 2) / (1 - w))
    obs = lad.loc["FG_NONISCH_CM"]
    diff = obs["b"] - b_non
    se_diff = float(np.sqrt(obs["se"] ** 2 + se_non ** 2))
    p_diff = float(2 * stats.norm.sf(abs(diff / se_diff)))

    rows = [
        dict(quantity="all heart failure (measured)", b=hf["b"], se=hf["se"],
             lo=hf["lo"], hi=hf["hi"], p=hf["p"], cases=CASES_ALL_HF, source="module 05"),
        dict(quantity="heart failure + coronary disease (measured)", b=chd["b"], se=chd["se"],
             lo=chd["lo"], hi=chd["hi"], p=chd["p"], cases=CASES_HF_CHD, source="module 05"),
        dict(quantity="non-ischaemic cardiomyopathy (measured)", b=obs["b"], se=obs["se"],
             lo=obs["lo"], hi=obs["hi"], p=obs["p"], cases=int(obs["cases"]),
             source="module 14"),
        dict(quantity="heart failure without coronary disease (implied by the mixture)",
             b=b_non, se=se_non, lo=b_non - 1.96 * se_non, hi=b_non + 1.96 * se_non,
             p=float(2 * stats.norm.sf(abs(b_non / se_non))),
             cases=CASES_ALL_HF - CASES_HF_CHD, source="arithmetic, w = 0.466"),
    ]
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False, float_format="%.6g")

    print(f"mixing weight w = {CASES_HF_CHD:,} / {CASES_ALL_HF:,} = {w:.3f}\n")
    for r in rows:
        print(f"  {r['quantity']:58s} {r['b']:+.4f} ({r['lo']:+.4f}, {r['hi']:+.4f}) "
              f"p={r['p']:.3g}")
    print(f"\n  implied minus observed: {diff:+.4f} (se {se_diff:.4f}), p={p_diff:.2f}")
    print("  -> " + ("consistent: the three numbers behave like one effect diluted by the "
                     "ischaemic half of heart failure" if p_diff > 0.05
                     else "NOT consistent: the mixture reading fails"))
    print(f"\n-> {OUT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

