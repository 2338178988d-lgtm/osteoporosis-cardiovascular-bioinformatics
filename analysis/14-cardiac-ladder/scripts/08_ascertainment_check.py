"""Is the Finnish excess a case-ascertainment artefact or a population effect?

Section H found the Finnish estimate is larger than the non-Finnish one in all three endpoints
measured in both settings (ratio 1.7, 2.2, 11.8; one significant). Three explanations fit
equally well so far:

  population       Finnish founder effects genuinely differ
  ascertainment    FinnGen endpoints are hospital ICD codes from a national registry, which
                   select more severe disease than the PheCode/UKB-coded comparators; a
                   severity-selected case mix inflates a per-SD effect
  winner's curse   FinnGen was analysed first in every endpoint

One cheap check separates a *general* Finnish inflation from endpoint-specific noise: apply the
same instruments to Finnish outcomes that are **not** expected to respond to bone density, and
see whether they are also inflated. If Finnish case ascertainment inflates everything, the
negative-control estimates should be systematically larger in Finland too. If only the
protective endpoints are inflated, the pattern is about the phenotype rather than the pipeline.

Negative controls already run in module 13: six FinnGen endpoints unrelated to bone density.
Their estimate magnitudes are tiny either way, so this check can only bound the mechanism, not
identify it - which is stated in the output rather than left implied.

Outputs (results/tables/14-cardiac-ladder/ascertainment_check.csv)

Usage:
  python analysis/14-cardiac-ladder/scripts/08_ascertainment_check.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / "results" / "tables" / "14-cardiac-ladder"


def main():
    neg = ROOT / "results" / "tables" / "13-robustness" / "negative_controls.csv"
    if not neg.exists():
        print("negative controls absent; run module 13 first")
        return 2
    d = pd.read_csv(neg)
    e = d[d["exposure"] == "eBMD"]
    # the module-13 table carries a heart-failure positive control alongside the six negative
    # controls; counting it as a null would put the very effect under study on both sides of
    # the comparison and inflate the null bound to exactly the protective estimate
    if "role" in e.columns:
        n_pos = int((e["role"] != "negative_control").sum())
        e = e[e["role"] == "negative_control"]
        if n_pos:
            print(f"  excluded {n_pos} non-negative-control row(s) "
                  f"(positive control) from the null side")
    if e.empty:
        print("no eBMD negative-control rows")
        return 2

    rows = []
    for _, s in e.iterrows():
        rows.append(dict(endpoint=s["label"], b=float(s["b"]), se=float(s["se"]),
                         p=float(s["p"]), power=float(s.get("power_ref", np.nan)),
                         role="negative control"))
    # the three endpoints whose Finnish estimate exceeds the non-Finnish one
    allc = pd.read_csv(OUT / "big_cohorts_all.csv")
    for endpoint, fin_label in [("AF", "FinnGen"), ("HF", "FinnGen"), ("CM", "FinnGen")]:
        d2 = allc[(allc["endpoint"] == endpoint)
                  & (allc["population"].str.startswith("Finnish"))]
        for _, s in d2.iterrows():
            rows.append(dict(endpoint=f"{endpoint}: {s['label']}", b=float(s["b"]),
                             se=float(s["se"]), p=np.nan, power=np.nan,
                             role="protective (Finnish)"))

    df = pd.DataFrame(rows)
    df["abs_b"] = df["b"].abs()
    df["abs_t"] = (df["b"] / df["se"]).abs()
    df.to_csv(OUT / "ascertainment_check.csv", index=False, float_format="%.6g")

    nc = df[df["role"] == "negative control"]
    pr = df[df["role"] == "protective (Finnish)"]
    print("=== Finnish estimates, negative controls vs protective endpoints")
    print(f"  negative controls (k={len(nc)}): mean |b| = {nc['abs_b'].mean():.4f}, "
          f"max |b| = {nc['abs_b'].max():.4f}, median |z| = {nc['abs_t'].median():.2f}")
    print(f"  protective       (k={len(pr)}): mean |b| = {pr['abs_b'].mean():.4f}, "
          f"min |b| = {pr['abs_b'].min():.4f}, median |z| = {pr['abs_t'].median():.2f}")
    overlap = pr["abs_b"].min() <= nc["abs_b"].max()
    print(f"\n  the smallest protective |b| is {pr['abs_b'].min():.4f}; "
          f"the largest negative-control |b| is {nc['abs_b'].max():.4f}")
    if overlap:
        print("  -> a general Finnish inflation cannot be excluded: some null endpoints move "
              "by as much as the smallest protective effect")
    else:
        print("  -> the protective endpoints are separated from every null endpoint; a general "
              "ascertainment or pipeline inflation is not visible at this precision")
    print("\n  limitation: negative controls are null by construction, so this check bounds the")
    print("  mechanism rather than identifying it. It cannot separate founder-population from")
    print("  case-severity selection within the protective endpoints themselves.")
    print(f"\n-> {(OUT / 'ascertainment_check.csv').relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

