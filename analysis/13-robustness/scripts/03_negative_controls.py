"""Module 13C 鈥?negative controls: does the instrument set have generic pleiotropy?

"Endpoint-specific" has so far been supported by four cardiovascular outcomes, which are
all vascular or cardiac and therefore not independent tests of specificity. The standard
check is a set of endpoints that bone density has no business influencing.

Six FinnGen R12 endpoints, pre-specified before the run and well powered (26,894 to 83,886
cases each, six organ systems):

  K11_APPENDACUT      acute appendicitis          (gastrointestinal / surgical)
  H7_CATARACTSENILE   senile cataract             (eye)
  L12_ATOPIC          atopic dermatitis           (skin)
  G6_MIGRAINE         migraine                    (neurological)
  J10_TONSILLECTOMY   tonsillectomy and adenotomy (ENT / practice pattern)
  J10_SINUSITIS       acute sinusitis             (ENT / infectious)

Decision rule, fixed before the run:
  0 of 6 significant at Bonferroni 0.05/6  -> the instrument set is not generically
                                              pleiotropic; "endpoint-specific" may be written
  1-2 nominal (p<0.05) but none Bonferroni  -> expected by chance for six tests (0.3
                                              expected); report the count and the threshold
  any Bonferroni-significant                -> generic pleiotropy; the main estimate needs a
                                              stricter instrument subset and the specificity
                                              claim must be withdrawn

Outputs (results/tables/13-robustness/):
  negative_controls.csv   one row per endpoint, with power and the Bonferroni threshold
  negative_summary.csv    count exceedances vs expectation, binomial test

Usage:
  python analysis/13-robustness/scripts/03_negative_controls.py
  python analysis/13-robustness/scripts/03_negative_controls.py --self-test
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from _common import estimate, have, mr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / "results" / "tables" / "13-robustness"
EXPOSURE = "eBMD"
EXPOSURES = ["eBMD", "FNBMD"]          # both, so the two BMD modalities can be compared
CONTROLS = {
    "FN_NEG_APPENDACUT": "acute appendicitis",
    "FN_NEG_CATARACT": "senile cataract",
    "FN_NEG_ATOPIC": "atopic dermatitis",
    "FN_NEG_MIGRAINE": "migraine",
    "FN_NEG_TONSILL": "tonsillectomy / adenotomy",
    "FN_NEG_SINUSITIS": "acute sinusitis",
}
# the two positive controls that must stay significant, for contrast in one table
POSITIVES = {"FG_HF": "heart failure (positive control)"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for exp in EXPOSURES:
        inst = mr.load_instruments(exp)
        if inst is None or inst.empty:
            print(f"{exp}: no instruments, skipped")
            continue
        print(f"\n=== exposure {exp}: {len(inst)} instruments")
        for key, label in {**POSITIVES, **CONTROLS}.items():
            r = estimate(exp, key, inst=inst, label=f"{exp} -> {label}")
            if r:
                r["role"] = "positive_control" if key in POSITIVES else "negative_control"
                r["label"] = label
                rows.append(r)
    if not rows:
        print("no endpoint could be read 鈥?are the files downloaded?")
        return 2

    df = pd.DataFrame(rows)
    summaries = []
    for exp, sub in df.groupby("exposure"):
        neg = sub[sub.role == "negative_control"]
        m = len(neg)
        thr = 0.05 / m if m else np.nan
        df.loc[neg.index, "bonferroni_threshold"] = thr
        n_nominal = int((neg["p"] < 0.05).sum())
        summaries.append(dict(
            exposure=exp, n_negative=m, n_nominal_p05=n_nominal,
            expected_nominal=0.05 * m, bonferroni_threshold=thr,
            n_passing_bonferroni=int(neg["p"].lt(thr).sum()),
            binomial_p_upper=float(stats.binom.sf(n_nominal - 1, m, 0.05)) if m else np.nan,
            min_p=float(neg["p"].min()) if m else np.nan,
            median_power=float(neg["power_ref"].median()) if m else np.nan))
        print(f"\n=== {exp}: {m} negative controls, Bonferroni threshold {thr:.4f}")
        for _, r in neg.sort_values("p").iterrows():
            flag = ("PASS*" if r["p"] < thr else "nominal" if r["p"] < 0.05 else "null")
            print(f"  {r['label']:30s} n={r['nsnp']:4d} b={r['b']:+.4f} p={r['p']:.3g} "
                  f"power={r['power_ref']:.0%}  {flag}")
    df["bonferroni_threshold"] = df["bonferroni_threshold"].where(
        df["bonferroni_threshold"].notna(),
        df["role"].map(lambda r: 0.05 if r == "positive_control" else np.nan))
    df["passes_bonferroni"] = df["p"] < df["bonferroni_threshold"]
    df.sort_values(["exposure", "role", "p"]).reset_index(drop=True).to_csv(
        OUT / "negative_controls.csv", index=False, float_format="%.6g")
    pd.DataFrame(summaries).to_csv(OUT / "negative_summary.csv", index=False,
                                   float_format="%.6g")
    for s in summaries:
        print(f"\n{s['exposure']}: nominal p<0.05 {s['n_nominal_p05']}/{s['n_negative']} "
              f"(expected {s['expected_nominal']:.2f}); binomial upper p="
              f"{s['binomial_p_upper']:.3g}; min p {s['min_p']:.3g}; "
              f"median power {s['median_power']:.0%}")
    print(f"\n-> {(OUT / 'negative_controls.csv').relative_to(ROOT).as_posix()}")
    return 0


def selftest():
    """Six null endpoints must yield no Bonferroni exceedance; a planted one must be caught.

    Scales match the real data (bx ~ 0.04, sx ~ 0.003, outcome SE ~ 0.02), so an IVW SE of
    ~0.02 makes a planted 0.08 detectable and a null not.
    """
    rng = np.random.default_rng(9)
    n = 592
    bx = rng.normal(0.04, 0.02, n)
    sx = np.full(n, 0.003)
    sy = np.full(n, 0.02)
    hits = 0
    for _ in range(6):
        hits += int(mr.ivw(bx, rng.normal(0, sy), sx, sy, random=True)["p"] < 0.05 / 6)
    planted = mr.ivw(bx, -0.08 * bx + rng.normal(0, sy), sx, sy, random=True)
    print(f"selftest: {hits}/6 null endpoints passed Bonferroni; planted b="
          f"{planted['b']:+.4f} p={planted['p']:.3g}")
    ok = hits == 0 and planted["p"] < 0.05 and planted["b"] < 0
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

