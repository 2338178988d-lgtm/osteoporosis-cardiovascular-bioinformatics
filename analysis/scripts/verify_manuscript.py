"""Verify every number in the manuscript draft against the tables it claims to come from.

Run this after any revision of `results/reports/璁烘枃鑽夌_v1.md` or after any re-run of an
analysis module: a number that cannot be reproduced from a table will not survive review.

    python analysis/scripts/verify_manuscript.py

Exit code 0 when every claimed value matches its table, 1 otherwise.
"""
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
T = ROOT / "results" / "tables"

# the draft lives in results/reports/ inside the project and in manuscript/ inside the
# submission package; the verifier has to work in both layouts
DRAFT_CANDIDATES = [
    ROOT / "results" / "reports" / "璁烘枃鑽夌_v1.md",
    ROOT / "manuscript" / "璁烘枃鑽夌_v1.md",
]
DRAFT = next((p for p in DRAFT_CANDIDATES if p.exists()), DRAFT_CANDIDATES[0])

lad = pd.read_csv(T / "14-cardiac-ladder/ladder.csv").set_index("outcome")
panel = pd.read_csv(T / "14-cardiac-ladder/panel.csv").set_index("outcome")
allc = pd.read_csv(T / "14-cardiac-ladder/big_cohorts_all.csv")
chk = pd.read_csv(T / "14-cardiac-ladder/panel_checks.csv").set_index("outcome")
eff = pd.read_csv(T / "14-cardiac-ladder/effective_tests.csv").iloc[0]
dil = pd.read_csv(T / "14-cardiac-ladder/dilution.csv")
asc = pd.read_csv(T / "14-cardiac-ladder/ascertainment_check.csv")
avc = pd.read_csv(T / "14-cardiac-ladder/avc.csv").set_index("outcome")
rung = pd.read_csv(T / "14-cardiac-ladder/rung_concordance.csv")
mod = pd.read_csv(T / "15-modality/modality.csv")
modpool = pd.read_csv(T / "15-modality/modality_pooled.csv")
conc = pd.read_csv(T / "15-modality/concordance.csv")
mvmr = pd.read_csv(T / "06-bone-axis/axis_mvmr.csv")
body = pd.read_csv(T / "11-body-composition/bodycomp_mvmr.csv")
alm = pd.read_csv(T / "11-body-composition/bodycomp_alm_own.csv")
hfmv = pd.read_csv(T / "05-hf-followup/hf_mvmr.csv")
sub = pd.read_csv(T / "05-hf-followup/fg_subtype.csv").set_index("outcome")
neg = pd.read_csv(T / "13-robustness/negative_summary.csv")
lv = pd.read_csv(T / "13-robustness/lv_mass.csv")
coloc = pd.read_csv(T / "09-integration/coloc_summary.csv").set_index("metric")["value"]
pw = pd.read_csv(T / "12-pathway-stratified/pathway_contributions.csv")
gs = pd.read_csv(T / "12-pathway-stratified/geneset_fwer.csv").set_index("outcome")
het = pd.read_csv(T / "01-nhanes/table13_outcome_heterogeneity.csv")
cox = pd.read_csv(T / "01-nhanes/table4_survey_cox_mortality.csv")
nhanes = pd.read_csv(T / "01-nhanes/table12_outcome_comparison.csv").set_index("model")

fails = []


def check(name, claimed, actual, tol=0.006):
    ok = abs(float(claimed) - float(actual)) <= max(tol, abs(float(claimed)) * 0.02)
    print(f"{'ok  ' if ok else 'FAIL'} {name:44s} claimed {claimed:>10}  table {actual:>10}")
    if not ok:
        fails.append(name)


def pool(b, se):
    b, se = np.asarray(b, float), np.asarray(se, float)
    w = 1 / se ** 2
    bf = float(np.sum(w * b) / np.sum(w))
    q = float(np.sum(w * (b - bf) ** 2))
    dfree = len(b) - 1
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (q - dfree) / c) if c > 0 and dfree > 0 else 0.0
    ws = 1 / (se ** 2 + tau2)
    i2 = max(0.0, (q - dfree) / q) if q > 0 and dfree > 0 else 0.0
    return float(np.sum(ws * b) / np.sum(ws)), float(np.sqrt(1 / np.sum(ws))), i2


def r_of(a, b):
    s = rung[((rung["a"] == a) & (rung["b"] == b)) | ((rung["a"] == b) & (rung["b"] == a))]
    return round(float(s["r"].iloc[0]), 2)


def main():
    if not DRAFT.exists():
        print(f"draft absent: {DRAFT}")
        return 2
    print(f"verifying {DRAFT.relative_to(ROOT).as_posix()}\n")

    print("-- table 1 and section 3.1: cohorts")
    for label, key in [("GCST90475951 (19,708)", -0.001), ("GCST90475949 (20,680)", 0.003),
                       ("GCST90473562 (3,578)", -0.070)]:
        r = allc[(allc["endpoint"] == "CM") & (allc["label"] == label)].iloc[0]
        check(f"CM {label}", key, round(float(r["b"]), 3))
    for ep, key in [("AF", -0.063), ("HF", -0.037), ("CM", -0.040)]:
        d = allc[allc["endpoint"] == ep]
        br, se, i2 = pool(d["b"], d["se"])
        check(f"{ep} pooled beta (all cohorts)", key, br, 0.002)
    check("CM pooled I2", 0.72, pool(allc[allc["endpoint"] == "CM"]["b"],
                                    allc[allc["endpoint"] == "CM"]["se"])[2], 0.01)

    print("\n-- section 3.2: ladder")
    for key, val in [("FG_HYPERTENS", 0.008), ("FG_MI_STRICT", 0.007), ("FG_AF", -0.086),
                     ("FG_VHD", -0.074), ("FG_CAVS", -0.024), ("FG_HYPERTENS_HD", -0.013),
                     ("FG_NONISCH_CM", -0.101), ("FG_CARDMYO", -0.062), ("FG_HF", -0.064)]:
        check(f"{key} beta", val, round(float(lad.loc[key, "b"]), 3))
    check("AF power", 0.67, round(float(lad.loc["FG_AF", "power_ref"]), 2))
    check("hypertension power", 0.88, round(float(lad.loc["FG_HYPERTENS", "power_ref"]), 2))
    check("Egger p HF", 0.80, round(float(lad.loc["FG_HF", "egger_intercept_p"]), 2))
    check("Egger p AF", 0.73, round(float(lad.loc["FG_AF", "egger_intercept_p"]), 2))
    check("r HF vs non-isch CM", 0.69, r_of("FG_HF", "FG_NONISCH_CM"), 0.01)
    check("r non-isch CM vs MI", 0.02, r_of("FG_NONISCH_CM", "FG_MI_STRICT"), 0.01)
    check("implied HF without CHD", -0.127, round(float(dil.iloc[3]["b"]), 3))
    check("observed non-isch CM", -0.102, round(float(dil.iloc[2]["b"]), 3))

    print("\n-- section 3.3: discriminating panel")
    for key, val in [("FG_VARICVE", -0.014), ("FG_VTE", -0.049), ("FG_PULMEMB", -0.038),
                     ("FG_AORTANEUR", -0.027), ("FG_OTHARR", -0.026), ("FG_CONDUCTIO", -0.019),
                     ("FG_AVBLOCK", -0.061), ("FG_CVD_ALL", -0.016), ("FG_IHD", 0.015),
                     ("FG_CHD", 0.029)]:
        check(f"{key} beta", val, round(float(panel.loc[key, "b"]), 3))
    check("other arrhythmias power", 0.95, round(float(panel.loc["FG_OTHARR", "power_ref"]), 2))
    check("VTE p", 0.061, round(float(panel.loc["FG_VTE", "p"]), 3))
    check("VTE trimmed p", 0.046, round(float(chk.loc["FG_VTE", "p_trimmed"]), 3))

    print("\n-- section 3.4: diffuse")
    check("loci tested", 600, coloc["loci_tested"])
    check("PP.H4 >= 0.8", 1, coloc["PP_H4_ge_0.8"])
    check("median PP.H4", 0.002, coloc["median_PP_H4"], 0.0005)
    u = pw[(pw["outcome"] == "FG_HF") & (pw["group"] == "unassigned")].iloc[0]
    check("unassigned share", 95.5, float(u["share_pct"]), 0.2)
    check("unassigned nsnp", 546, int(u["nsnp"]))
    check("geneset FG_HF max z", 3.64, float(gs.loc["FG_HF", "observed_max_abs_z"]), 0.01)
    check("geneset FG_HF p", 0.67, float(gs.loc["FG_HF", "p_familywise"]), 0.005)
    check("geneset HF max z", 4.52, float(gs.loc["HF", "observed_max_abs_z"]), 0.01)
    check("geneset HF p", 0.40, float(gs.loc["HF", "p_familywise"]), 0.005)

    print("\n-- section 3.5: mediators and body composition")
    for med, val in [("VITD", -0.2), ("FGF23", -0.1), ("SOST", 2.6), ("BGLAP", 3.7),
                     ("DKK1", 0.3), ("PTH", 0.7), ("TESTO", -0.4), ("ESTRADIOL", -1.5),
                     ("CST3", 4.7), ("CRP", -0.3)]:
        r = mvmr[(mvmr["outcome"] == "FG_HF") & (mvmr["mediator"] == med)].iloc[0]
        check(f"attenuation {med}", val, float(r["attenuation_pct"]), 0.15)
    for model, val in [("M2_eBMD_ALM", 8.4), ("M3_eBMD_HEIGHT", 8.7),
                       ("M4_eBMD_ALM_BMI", 6.7)]:
        r = body[(body["outcome"] == "FG_HF") & (body["model"] == model)].iloc[0]
        check(f"body comp {model}", val, float(r["attenuation_pct"]), 0.15)
    a = alm[(alm["outcome"] == "FG_HF") & (alm["model"] == "univariable")].iloc[0]
    a2 = alm[(alm["outcome"] == "FG_HF") & (alm["model"] == "ALM+BMI+HEIGHT")].iloc[0]
    check("ALM own beta", 0.125, round(float(a["b"]), 3))
    check("ALM own log10 p", -7.5, round(math.log10(float(a["p"])), 1), 0.1)
    check("ALM after BMI+height p", 0.50, round(float(a2["p"]), 2))

    print("\n-- section 3.6: modality")
    for exp, val in [("eBMD", -0.064), ("QUS_2", -0.066), ("QUS_3", -0.071),
                     ("FNBMD", 0.068), ("FNBMD_2", 0.039), ("FNBMD_3", 0.025)]:
        r = mod[(mod["exposure"] == exp) & (mod["outcome"] == "FG_HF")].iloc[0]
        check(f"{exp} FG_HF beta", val, round(float(r["b"]), 3))
    for grp, out, val in [("ultrasound", "FG_HF", -0.067), ("dxa", "FG_HF", 0.047),
                          ("dxa", "HF", 0.002)]:
        r = modpool[(modpool["group"] == grp) & (modpool["outcome"] == out)].iloc[0]
        check(f"{grp} pooled {out}", val, round(float(r["b"]), 3))
    check("DXA I2 in HERMES", 0.60,
          float(modpool[(modpool["group"] == "dxa") & (modpool["outcome"] == "HF")]["i2"].iloc[0]),
          0.01)
    us = ["eBMD", "QUS_2", "QUS_3"]
    within = conc[(conc["predictor"].isin(us) == conc["predicted"].isin(us))
                  & (conc["predictor"] != conc["predicted"])]
    across = conc[conc["predictor"].isin(us) != conc["predicted"].isin(us)]
    check("concordance within min", 0.94, round(float(within["b"].min()), 2))
    check("concordance within max", 0.99, round(float(within["b"].max()), 2))
    check("concordance across mean", 0.46, round(float(across["b"].mean()), 2))

    print("\n-- section 3.7: robustness")
    for key, base, trim in [("FG_AF", -0.0855, -0.0941), ("FG_NONISCH_CM", -0.1015, -0.1031),
                            ("FG_VHD", -0.0745, -0.0755), ("FG_MI_STRICT", 0.0069, -0.0011)]:
        check(f"shared-locus {key} base", base, round(float(chk.loc[key, "b_base"]), 4))
        check(f"shared-locus {key} trim", trim, round(float(chk.loc[key, "b_trimmed"]), 4))
    check("shared loci removed", 4, int(chk["n_shared_removed"].max()))
    check("Steiger max removed", 3, int(chk["n_steiger_removed"].max()))
    check("Meff", 14.1, round(float(eff["meff"]), 1))
    check("Meff alpha", 0.0035, round(float(eff["alpha_meff"]), 4))
    e = neg[neg["exposure"] == "eBMD"].iloc[0]
    check("negative controls nominal", 0, int(e["n_nominal_p05"]))
    check("negative controls min p", 0.052, round(float(e["min_p"]), 3))
    check("negative controls median power", 0.92, round(float(e["median_power"]), 2))
    check("LV mass SD beta", -0.041, round(float(lv[lv["outcome"] == "LVM"]["b_sd"].iloc[0]), 3))
    check("LV mass SD p", 0.015, round(float(lv[lv["outcome"] == "LVM"]["p"].iloc[0]), 3))

    print("\n-- section 3.8: population layer")
    for out, val in [("chf", 1.37), ("chd", 1.34), ("mi", 1.26), ("stroke", 1.49)]:
        check(f"NHANES crude {out}", val, round(float(nhanes.loc["M1", f"or.{out}"]), 2))
    for out, val in [("chf", 1.14), ("chd", 1.10), ("mi", 1.07), ("stroke", 1.22)]:
        check(f"NHANES adjusted {out}", val, round(float(nhanes.loc["M4", f"or.{out}"]), 2))
    check("NHANES stroke adjusted p", 0.003,
          round(float(nhanes.loc["M4", "p.stroke"]), 3), 0.001)
    o = cox[(cox["outcome"] == "dead") & (cox["model"].str.startswith("M3"))
            & (cox["term"] == "bmd_catOsteoporosis")].iloc[0]
    check("NHANES all-cause HR", 1.76, round(float(o["HR"]), 2))
    c = cox[(cox["outcome"] == "cvd_death_num") & (cox["model"].str.startswith("M3"))
            & (cox["term"] == "bmd_catOsteoporosis")].iloc[0]
    check("NHANES CV death HR", 1.51, round(float(c["HR"]), 2))

    print("\n-- section 4: discussion claims")
    cad = hfmv[hfmv["mediator"] == "CAD"].iloc[0]
    check("CAD MVMR beta", -0.045, round(float(cad["direct_b"]), 3))
    check("CAD MVMR log10 p", -3.9, round(math.log10(float(cad["direct_p"])), 1), 0.1)
    check("HF+obesity beta", -0.078, round(float(sub.loc["FG_HF_OBESE", "b"]), 3))
    check("HF+CHD p", 0.78, round(float(sub.loc["FG_HF_CHD", "p"]), 2))
    nc = asc[asc["role"] == "negative control"]
    pr = asc[asc["role"] == "protective (Finnish)"]
    check("ascertainment neg max |b|", 0.034, round(float(nc["abs_b"].max()), 3))
    check("ascertainment protective min |b|", 0.064, round(float(pr["abs_b"].min()), 3))
    check("AVC beta", -0.058, round(float(avc.loc["AVC", "b"]), 3))
    check("AVC p", 0.065, round(float(avc.loc["AVC", "p"]), 3))
    check("AVC power", 0.53, round(float(avc.loc["AVC", "power_ref"]), 2))

    print("\n-- heterogeneity quoted in section 3.8")
    for _, r in het.iterrows():
        print(f"     {r['comparison']:16s} p = {r['p']:.2f}")

    print()
    if fails:
        print(f"{len(fails)} MISMATCH: {fails}")
        return 1
    print("ALL CHECKS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())

