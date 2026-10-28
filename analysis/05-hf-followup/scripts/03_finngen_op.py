"""FinnGen osteoporosis as an independent exposure: does OP affect HF, and does HF
affect OP, when both sides are external to UK Biobank?

The UK-Biobank OP instrument set has only 19 variants and collapses under any
sensitivity analysis. FinnGen R12 `M13_OSTEOPOROSIS` (10,461 cases / 473,264 controls)
is a different population, so it provides an OP instrument set that shares no sample
with eBMD.

Outputs (results/tables/05-hf-followup/fg_op.csv):
  OP -> heart failure, both inside FinnGen (same population, case-control overlap noted)
  OP -> eBMD (UK Biobank), a clean cross-population test
  heart failure -> OP, within FinnGen

Usage: python analysis/05-hf-followup/scripts/03_finngen_op.py
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
_spec = importlib.util.spec_from_file_location("fg", HERE / "02_finngen.py")
fg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fg)
mr = fg.mr

OUT = fg.OUT


def run(exp_key, out_key, inst=None, label=""):
    if inst is None:
        inst = mr.load_instruments(exp_key)
    if inst is None or inst.empty:
        print(f"  {label or exp_key}: no instruments")
        return None
    out = mr.outcome_at_snps(out_key, inst["SNP"].tolist())
    if out is None or out.empty:
        print(f"  {label or exp_key}: outcome has no matching SNPs")
        return None
    d = mr.harmonise(inst, out)
    if len(d) < 5:
        print(f"  {label or exp_key}: only {len(d)} usable SNPs")
        return None
    bx, by = d["bx"].to_numpy(float), d["by"].to_numpy(float)
    sx, sy = d["sx"].to_numpy(float), d["sy"].to_numpy(float)
    r = mr.ivw(bx, by, sx, sy, random=True)
    e = mr.egger(bx, by, sx, sy) if len(d) >= 10 else dict(intercept_p=np.nan, b=np.nan, p=np.nan)
    row = dict(pair=label or f"{exp_key} -> {out_key}", exposure=exp_key, outcome=out_key,
               nsnp=len(d), b=r["b"], se=r["se"], p=r["p"],
               lo=r["b"] - 1.96 * r["se"], hi=r["b"] + 1.96 * r["se"],
               or_=float(np.exp(r["b"])),
               egger_intercept_p=e["intercept_p"], q_p=r["q_p"],
               mean_F=float((d["bx"] / d["sx"]).pow(2).mean()),
               sum_r2=float((d["bx"].pow(2) / (d["bx"].pow(2) + d["sx"].pow(2) * d["n_x"])).sum()))
    print(f"  {row['pair']:34s} n={len(d):4d}  b={r['b']:+.4f} "
          f"({row['lo']:+.4f}, {row['hi']:+.4f})  p={r['p']:.3g}  OR={row['or_']:.3f}  "
          f"Egger int p={e['intercept_p']:.3g}")
    return row


def main():
    print("== FinnGen osteoporosis instrument set")
    fg_inst = fg.build_finngen_instruments("FG_OP")
    if fg_inst is None or fg_inst.empty:
        print("  could not build FG_OP instruments")
        return 1

    rows = []
    print("\n== OP -> HF in FinnGen (same population as the exposure GWAS)")
    for key, lab in [("FG_HF", "FG-OP -> FG-HF (strict)"),
                     ("FG_HF_CHD", "FG-OP -> FG-HF + CHD"),
                     ("FG_HF_OBESE", "FG-OP -> FG-HF + BMI>=25")]:
        r = run("FG_OP", key, inst=fg_inst, label=lab)
        if r:
            rows.append(r)

    print("\n== OP -> eBMD across populations (FinnGen exposure, UKB outcome)")
    r = run("FG_OP", "eBMD", inst=fg_inst, label="FG-OP -> UKB eBMD")
    if r:
        rows.append(r)

    print("\n== HF -> OP within FinnGen (exposure instruments from FG_HF)")
    hf_inst = fg.build_finngen_instruments("FG_HF")
    r = run("FG_HF", "FG_OP", inst=hf_inst, label="FG-HF -> FG-OP")
    if r:
        rows.append(r)

    if rows:
        pd.DataFrame(rows).to_csv(OUT / "fg_op.csv", index=False, float_format="%.6g")
        print(f"\n-> {OUT / 'fg_op.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

