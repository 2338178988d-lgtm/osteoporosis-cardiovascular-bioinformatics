"""Module 06 鈥?does the bone-to-heart-failure effect travel through a bone-derived factor?

The eBMD -> heart-failure effect survives every test thrown at it (module 05), but the
pathway is unknown. The mediators tested so far (BMI, type-2 diabetes, CAD) are *not*
bone-derived, so ruling them out only says the effect is not routed through metabolic or
atherosclerotic disease. This module tests the factors a bone cell actually secretes and
that are known to act on the heart:

  FGF23     osteocyte-derived; drives cardiomyocyte hypertrophy via FGFR4
  SOST      sclerostin; Wnt antagonist, bone formation and cardiac fibrosis
  BGLAP     osteocalcin; undercarboxylated form acts on myocardial energy metabolism
  CST3      cystatin C; renal-function proxy AND heart-failure prognostic marker
  CRP       inflammation; the leading competing mechanism for heart failure

Design: multivariable MR of eBMD and one mediator at a time on the heart-failure
outcome, using eBMD's LD-clumped instruments. Three quantities are reported per
mediator:

  direct effect    eBMD effect on HF adjusting for the mediator (MVMR-IVW, sandwich SE)
  attenuation      how much of the univariable eBMD effect the mediator absorbs
  mediator effect  the mediator's own conditional effect on HF

Plus, for context, univariable MR of each mediator on HF using that mediator's own
instruments, which shows whether the mediator moves HF at all.

Outputs (results/tables/06-bone-axis/):
  axis_mvmr.csv        eBMD direct effect and attenuation for each mediator
  axis_univariable.csv mediator -> HF on the mediator's own instruments
  axis_strength.csv    instrument strength and conditional F for each mediator

Usage:
  python analysis/06-bone-axis/scripts/01_axis_mvmr.py
  python analysis/06-bone-axis/scripts/01_axis_mvmr.py --self-test
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
_fu = importlib.util.spec_from_file_location("fu", HERE.parents[1] / "05-hf-followup" / "scripts" / "01_hf_followup.py")
fu = importlib.util.module_from_spec(_fu)
_fu.loader.exec_module(fu)

ROOT = HERE.parents[2]
OUT = ROOT / "results" / "tables" / "06-bone-axis"

EXPOSURE = "eBMD"
# two heart-failure outcomes: the overlapping discovery set and the independent one
OUTCOMES = {"HF": "HERMES (overlapping)", "FG_HF": "FinnGen R12 (independent)"}
MEDIATORS = ["FGF23", "SOST", "BGLAP", "CST3", "CRP", "DKK1", "PTH", "TESTO", "ESTRADIOL", "VITD"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()

    OUT.mkdir(parents=True, exist_ok=True)
    mv_rows, uv_rows, st_rows = [], [], []

    for out_key, out_label in OUTCOMES.items():
        print(f"\n=== outcome: {out_key} ({out_label})")
        tab = fu.stack(EXPOSURE, out_key, MEDIATORS)
        if tab is None:
            print("  no data")
            continue
        present = [m for m in MEDIATORS if f"b_{m}" in tab.columns and tab[f"b_{m}"].notna().any()]
        print(f"  {len(tab)} instruments; mediators with usable effects: {', '.join(present)}")

        base_bx = tab["bx"].to_numpy(float)
        base_by = tab["by"].to_numpy(float)
        base_sx = tab["sx"].to_numpy(float)
        base_sy = tab["sy"].to_numpy(float)
        base = mr.ivw(base_bx, base_by, base_sx, base_sy, random=True)
        print(f"  univariable eBMD -> {out_key}: b={base['b']:+.4f} (p={base['p']:.3g})")

        for med in present:
            sub = tab.dropna(subset=[f"b_{med}", f"s_{med}"])
            if len(sub) < 20:
                print(f"  {med}: only {len(sub)} shared instruments, skipped")
                continue
            X = np.column_stack([sub["bx"].to_numpy(float), sub[f"b_{med}"].to_numpy(float)])
            y = sub["by"].to_numpy(float)
            sy = sub["sy"].to_numpy(float)
            res = fu.mvmr_ivw(X, y, sy)
            if res is None:
                print(f"  {med}: singular design")
                continue
            beta, se, resid = res
            p = 2 * stats.norm.sf(np.abs(beta / se))
            cf = fu.conditional_f(sub["bx"].to_numpy(float), sub["sx"].to_numpy(float),
                                  sub[f"b_{med}"].to_numpy(float), len(sub))
            corr = float(np.corrcoef(sub["bx"].to_numpy(float), sub[f"b_{med}"].to_numpy(float))[0, 1])
            # the attenuation baseline must come from the SAME instruments, otherwise the
            # comparison mixes a precise full-set estimate with a small subset estimate
            sub_base = mr.ivw(sub["bx"].to_numpy(float), y, sub["sx"].to_numpy(float), sy,
                              random=True)
            mv_rows.append(dict(outcome=out_key, mediator=med, nsnp=len(sub),
                                univariable_b=sub_base["b"], univariable_p=sub_base["p"],
                                univariable_full_b=base["b"], univariable_full_p=base["p"],
                                direct_b=beta[0], direct_se=se[0], direct_p=p[0],
                                direct_or=float(np.exp(beta[0])),
                                direct_or_lo=float(np.exp(beta[0] - 1.96 * se[0])),
                                direct_or_hi=float(np.exp(beta[0] + 1.96 * se[0])),
                                mediator_b=beta[1], mediator_se=se[1], mediator_p=p[1],
                                attenuation_pct=100 * (1 - beta[0] / sub_base["b"])))
            st_rows.append(dict(outcome=out_key, mediator=med, nsnp=len(sub),
                                conditional_F_eBMD=cf, corr_eBMD_mediator=corr,
                                resid_scale=float(np.sqrt(np.mean(resid ** 2)))))
            print(f"  {med:6s} n={len(sub):4d}  eBMD direct {beta[0]:+.4f} "
                  f"(p={p[0]:.3g}) vs same-instrument univariable {sub_base['b']:+.4f} "
                  f"(p={sub_base['p']:.3g}), attenuation "
                  f"{100 * (1 - beta[0] / sub_base['b']):+6.1f}%; "
                  f"{med} conditional {beta[1]:+.4f} (p={p[1]:.3g})")

        # univariable mediator -> HF on the mediator's own instruments
        for med in present:
            inst = mr.load_instruments(med)
            if inst is None or inst.empty:
                continue
            out = mr.outcome_at_snps(out_key, inst["SNP"].tolist())
            if out is None or out.empty:
                continue
            d = mr.harmonise(inst, out)
            if len(d) < 10:
                print(f"  {med} -> {out_key}: only {len(d)} instruments, skipped")
                continue
            bx, by = d["bx"].to_numpy(float), d["by"].to_numpy(float)
            sx, sy = d["sx"].to_numpy(float), d["sy"].to_numpy(float)
            r = mr.ivw(bx, by, sx, sy, random=True)
            e = mr.egger(bx, by, sx, sy) if len(d) >= 10 else dict(intercept_p=np.nan, p=np.nan)
            uv_rows.append(dict(mediator=med, outcome=out_key, nsnp=len(d),
                                b=r["b"], se=r["se"], p=r["p"],
                                lo=r["b"] - 1.96 * r["se"], hi=r["b"] + 1.96 * r["se"],
                                or_=float(np.exp(r["b"])),
                                q_p=r["q_p"], egger_intercept_p=e["intercept_p"]))
            print(f"  {med} -> {out_key} (own instruments): n={len(d):4d} b={r['b']:+.4f} "
                  f"(p={r['p']:.3g})")

    if mv_rows:
        pd.DataFrame(mv_rows).to_csv(OUT / "axis_mvmr.csv", index=False, float_format="%.6g")
        print(f"\n-> {OUT / 'axis_mvmr.csv'}")
    if uv_rows:
        pd.DataFrame(uv_rows).to_csv(OUT / "axis_univariable.csv", index=False, float_format="%.6g")
        print(f"-> {OUT / 'axis_univariable.csv'}")
    if st_rows:
        pd.DataFrame(st_rows).to_csv(OUT / "axis_strength.csv", index=False, float_format="%.6g")
        print(f"-> {OUT / 'axis_strength.csv'}")
    return 0


def self_test():
    """A mediator that carries half the effect must show ~50 % attenuation; a null
    mediator must show ~0 %. Both are simulated with a known causal structure."""
    rng = np.random.default_rng(3)
    n = 400
    bx = rng.normal(0, 0.05, n)
    sx = np.full(n, 0.004)
    bm = 0.6 * bx + rng.normal(0, 0.05, n)      # mediator partly driven by BMD
    sm = np.full(n, 0.006)
    sy = np.full(n, 0.02)
    rng_null = np.random.default_rng(4)
    bm_null = rng_null.normal(0, 0.05, n)       # mediator independent of BMD
    sy_null = sy

    y_carried = 0.10 * bx + 0.50 * bm + rng.normal(0, sy)     # 0.10 direct, 0.30 via mediator
    y_null = 0.40 * bx + 0.50 * bm_null + rng.normal(0, sy_null)

    out = []
    for label, bm_use, y in [("carried", bm, y_carried), ("null", bm_null, y_null)]:
        base = mr.ivw(bx, y, sx, sy, random=True)
        X = np.column_stack([bx, bm_use])
        beta, se, _ = fu.mvmr_ivw(X, y, sy)
        att = 100 * (1 - beta[0] / base["b"])
        out.append((label, base["b"], beta[0], att))
        print(f"selftest {label}: univariable {base['b']:+.3f} -> direct {beta[0]:+.3f}, "
              f"attenuation {att:+.1f}%")
    carried_att = out[0][3]
    null_att = out[1][3]
    ok = carried_att > 40 and abs(null_att) < 25 and abs(out[0][2] - 0.10) < 0.12
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

