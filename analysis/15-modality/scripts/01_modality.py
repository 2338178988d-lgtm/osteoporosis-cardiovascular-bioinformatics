"""Module 15 鈥?is the heel-ultrasound effect an exposure artefact, or a modality effect?

Module 13 found that a DXA femoral-neck BMD GWAS gives an estimate of the **opposite sign**
to the heel-ultrasound eBMD estimate (+0.068 vs -0.064, FinnGen, non-overlapping CIs) and
that no diagnostic removed the disagreement. With one file on each side, two readings stayed
open and could not be separated:

    biological  quantitative ultrasound reflects trabecular architecture and elasticity,
                which areal BMD does not capture;
    technical   the DXA file is a secondary deposit from a 2024 machine-learning methods
                paper, and no independent DXA GWAS existed to check it against.

That ambiguity is resolvable with published data, because the confound is two-dimensional:
**modality** (ultrasound vs DXA) and **source** (which paper, which cohort, which pipeline).
One file per cell cannot tell them apart; several can.

    heel ultrasound      eBMD    GCST90029004  PMID 29892013  n=445,921   (current exposure)
                         QUS_2   GCST90025982  PMID 34226706  n=445,855   independent paper
                         QUS_3   GCST006979    PMID 30598549  n=426,824   independent paper
    DXA femoral neck     FNBMD   GCST90446641  PMID 39349818  n=79,940    methods paper
                         FNBMD_2 GCST90446628  PMID 39349818  n=47,205    same paper (replicate)
                         FNBMD_3 GCST90568458  PMID 41206123  n=32,017    independent paper

DECISION RULES, written before the run:

  1. every ultrasound file is protective and every DXA file is positive, with pooled
     within-modality I^2 low
        -> the effect is a property of the phenotype, not of one file. The wording
           "heel-ultrasound bone density" survives, and the modality contrast becomes a
           finding rather than a limitation.
  2. the files disagree *within* a modality
        -> source dominates modality; the disagreement is about the data, not the phenotype,
           and the exposure-specific wording must be withdrawn in favour of "unresolved".
  3. the DXA estimates are all positive but the pooled CI includes zero
        -> the reversal is not established; fall back to "DXA did not show protection" and
           demote the modality claim to a limitation.

Each exposure is also predicted by the others' instruments (reciprocal prediction). Identical
files must predict each other strongly and with the right sign; if a file's instruments fail
to predict the same phenotype measured elsewhere, that file is suspect regardless of what it
estimates for heart failure.

Outputs (results/tables/15-modality/):
  modality.csv     six exposures x two outcomes, with power and diagnostics
  concordance.csv  reciprocal prediction between exposures

Usage:
  python analysis/15-modality/scripts/01_modality.py
  python analysis/15-modality/scripts/01_modality.py --self-test
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
import _summary as summary  # noqa: E402  (module-15 local helper)

OUT = ROOT / "results" / "tables" / "15-modality"

# (key, modality, source, pmid, n, role)
EXPOSURES = [
    ("eBMD",     "ultrasound", "Morris 2018",     "29892013", 445_921, "current exposure"),
    ("QUS_2",    "ultrasound", "PMID 34226706",   "34226706", 445_855, "independent paper"),
    ("QUS_3",    "ultrasound", "PMID 30598549",   "30598549", 426_824, "independent paper"),
    ("FNBMD",    "dxa",        "PMID 39349818",   "39349818", 79_940,  "methods paper"),
    ("FNBMD_2",  "dxa",        "PMID 39349818",   "39349818", 47_205,  "same paper (replicate)"),
    ("FNBMD_3",  "dxa",        "PMID 41206123",   "41206123", 32_017,  "independent paper"),
]
OUTCOMES = [("FG_HF", "FinnGen R12 heart failure (no UKB overlap)"),
            ("HF", "HERMES heart failure (shares UKB)")]
# reciprocal prediction pairs: (exposure, predicted trait)
CONCORDANCE = [("eBMD", "QUS_2"), ("eBMD", "QUS_3"), ("QUS_2", "QUS_3"),
               ("FNBMD", "FNBMD_2"), ("FNBMD", "FNBMD_3"), ("FNBMD_2", "FNBMD_3"),
               ("eBMD", "FNBMD"), ("eBMD", "FNBMD_3"), ("QUS_2", "FNBMD_3")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for key, modality, source, pmid, n, role in EXPOSURES:
        inst = mr.load_instruments(key)
        if inst is None or inst.empty:
            print(f"{key}: no instruments 鈥?run 02_extract_instruments.py / 02b_ld_clump.py")
            continue
        print(f"=== {key} ({modality}, {source}): {len(inst)} instruments, "
              f"mean R2={inst['r2'].mean():.4g}, min F={inst['F'].min():.1f}")
        for out_key, label in OUTCOMES:
            r = estimate(key, out_key, inst=inst, label=f"{key} -> {out_key}")
            if r:
                r.update(modality=modality, source=source, pmid=pmid, n_exposure=n, role=role)
                rows.append(r)
    if not rows:
        print("nothing estimated")
        return 1
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "modality.csv", index=False, float_format="%.6g")
    print(f"\n-> {(OUT / 'modality.csv').relative_to(ROOT).as_posix()}")

    conc = []
    for exp, target in CONCORDANCE:
        inst = mr.load_instruments(exp)
        if inst is None or not have(target):
            print(f"  concordance {exp} -> {target}: input missing, skipped")
            continue
        r = estimate(exp, target, inst=inst, label=f"{exp} -> {target}")
        if r:
            r.update(pair_label=f"{exp}->{target}", predicted=target, predictor=exp)
            conc.append(r)
    if conc:
        pd.DataFrame(conc).to_csv(OUT / "concordance.csv", index=False, float_format="%.6g")
        print(f"-> {(OUT / 'concordance.csv').relative_to(ROOT).as_posix()}")
        summary.report_concordance(pd.DataFrame(conc))

    pools = []
    for out_key, _ in OUTCOMES:
        pools.append(summary.report(df, out_key))
    pooled = pd.concat(pools, ignore_index=True)
    pooled.to_csv(OUT / "modality_pooled.csv", index=False, float_format="%.6g")
    print(f"-> {(OUT / 'modality_pooled.csv').relative_to(ROOT).as_posix()}")
    return 0


def selftest():
    """Pooling must recover a planted per-file effect and report heterogeneity when present."""
    rng = np.random.default_rng(7)
    b = rng.normal(-0.06, 0.005, 3)
    se = np.full(3, 0.02)
    p = summary.pool(b, se)
    print(f"selftest pool b={p['b']:+.4f} se={p['se']:.4f} I2={p['i2']:.2f}")
    ok = abs(p["b"] + 0.06) < 0.01
    b2 = np.array([-0.06, +0.06, +0.06])
    p2 = summary.pool(b2, se)
    print(f"selftest heterogeneous I2={p2['i2']:.2f} Qp={p2['q_p']:.3g}")
    ok = ok and p2["i2"] > 0.5
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

