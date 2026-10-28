"""Module 13A 鈥?does the heart-failure effect belong to bone density, or to heel ultrasound?

Every genetic result in this project uses one exposure: UK Biobank heel ultrasound eBMD
(GCST90029004, n=445,921). Ultrasound speed-of-sound and broadband attenuation are not DXA
bone mineral content; they depend on the soft-tissue path and on heel geometry. A distinct
modality answers a wording question that no amount of re-analysis of the same file can:

    DXA femoral neck BMD (GCST90446641, n=79,940 European) -> heart failure

Estimator: the same IVW / MR-Egger used everywhere else, on instruments built by the same
pipeline (p<5e-8, MAF 1-99 %, PLINK clump r2<0.001 within 10 Mb, 1000G EUR).

Decision rule, fixed before the run:
  direction agrees and the CI excludes zero
        -> the effect belongs to bone density; heel ultrasound is only the measurement
  null with the CI excluding -0.064
        -> the finding is specific to the heel-ultrasound phenotype and must be named that
           way, with soft tissue as a candidate explanation
  null with the CI containing -0.064
        -> underpowered; the modality question stays open, report the power

Outputs (results/tables/13-robustness/):
  dxa_exposure.csv     the two estimates with power
  dxa_concordance.csv  agreement between eBMD and DXA instruments at shared SNPs

Usage:
  python analysis/13-robustness/scripts/01_dxa_exposure.py
  python analysis/13-robustness/scripts/01_dxa_exposure.py --self-test
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from _common import estimate, have, mr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / "results" / "tables" / "13-robustness"
EXPOSURE = "FNBMD"
REFERENCE_EXPOSURE = "eBMD"
OUTCOMES = {"FG_HF": "FinnGen R12 (no UKB overlap)",
            "HF": "HERMES (shares UK Biobank)"}


def concordance():
    """Compare the two exposures' instrument sets at the SNPs they share."""
    inst = mr.load_instruments(EXPOSURE)
    ref = mr.load_instruments(REFERENCE_EXPOSURE)
    if inst is None or ref is None:
        return None
    a = ref.set_index("SNP")[["effect_allele", "other_allele", "bx"]].rename(
        columns={"bx": "bx_ref"})
    b = inst.set_index("SNP")[["effect_allele", "other_allele", "bx"]].rename(
        columns={"bx": "bx_new"})
    joined = a.join(b, how="inner", lsuffix="_ref", rsuffix="_new")
    flip = ((joined["effect_allele_ref"] == joined["other_allele_new"])
            & (joined["other_allele_ref"] == joined["effect_allele_new"]))
    joined.loc[flip, "bx_new"] = -joined.loc[flip, "bx_new"]
    ok = (joined["effect_allele_ref"] == joined["effect_allele_new"]) | flip
    joined = joined[ok]
    if len(joined) < 5:
        return None
    r = float(np.corrcoef(joined["bx_ref"], joined["bx_new"])[0, 1])
    slope = float(np.polyfit(joined["bx_ref"], joined["bx_new"], 1)[0])
    row = dict(exposure=EXPOSURE, reference=REFERENCE_EXPOSURE, n_shared=len(joined),
               n_flipped=int(flip.sum()), corr_effects=r, slope=bx_slope(slope),
               sign_agreement=float(np.mean(np.sign(joined["bx_ref"])
                                            == np.sign(joined["bx_new"]))))
    print(f"  shared instruments: {len(joined)} ({int(flip.sum())} allele-flipped); "
          f"effect correlation r={r:.3f}, sign agreement {row['sign_agreement']:.0%}")
    return pd.DataFrame([row])


def bx_slope(s):
    return float(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    if not have(EXPOSURE):
        print(f"{EXPOSURE}: file absent 鈥?run the download first")
        return 2
    inst = mr.load_instruments(EXPOSURE)
    if inst is None or inst.empty:
        print(f"{EXPOSURE}: no instrument table 鈥?run 02_extract_instruments.py")
        return 2
    n, mean_r2, min_f = len(inst), float(inst["r2"].mean()), float(inst["F"].min())
    print(f"=== DXA femoral neck BMD: {n} instruments, mean R2={mean_r2:.3g}, "
          f"min F={min_f:.1f}")

    rows = []
    for key, label in OUTCOMES.items():
        r = estimate(EXPOSURE, key, inst=inst,
                     label=f"DXA femoral neck BMD -> {key} ({label})")
        if r:
            rows.append(r)
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(OUT / "dxa_exposure.csv", index=False, float_format="%.6g")
        print(f"\n-> {(OUT / 'dxa_exposure.csv').relative_to(ROOT).as_posix()}")

    conc = concordance()
    if conc is not None:
        conc.insert(0, "nsnp_new", n)
        conc.to_csv(OUT / "dxa_concordance.csv", index=False, float_format="%.6g")
        print(f"-> {(OUT / 'dxa_concordance.csv').relative_to(ROOT).as_posix()}")
    return 0


def selftest():
    """The estimator must recover a planted effect and leave a null outcome alone.

    This checks estimator correctness, not realistic power: the synthetic set is deliberately
    well powered (large per-SNP effects, small outcome SE) so that a failure means the code
    is wrong rather than that the simulation was underpowered. Realistic power for the DXA
    design is reported at the end of `main`, not here.
    """
    rng = np.random.default_rng(4)
    n = 200
    bx = rng.normal(0.06, 0.02, n)
    sx = np.full(n, 0.004)
    sy = np.full(n, 0.01)
    out = []
    for by, tag in [(-0.064 * bx + rng.normal(0, sy), "planted -0.064"),
                    (rng.normal(0, sy), "null")]:
        r = mr.ivw(bx, by, sx, sy, random=True)
        out.append((tag, r["b"], r["p"]))
        print(f"selftest {tag:16s} b={r['b']:+.4f} (se {r['se']:.4f}) p={r['p']:.3g}")
    b_sig, p_sig = out[0][1], out[0][2]
    p_null = out[1][2]
    ok = abs(b_sig + 0.064) < 0.02 and p_sig < 0.01 and p_null > 0.01
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

