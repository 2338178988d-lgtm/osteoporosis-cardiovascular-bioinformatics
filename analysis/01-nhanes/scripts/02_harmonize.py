"""Step 2 鈥?harmonise NHANES 2005-2018 components into one analysable table.

Study question: is low bone mineral density (BMD) associated with prevalent
cardiovascular disease (CVD), and with subsequent all-cause / cardiovascular
mortality, in US adults aged >= 40?

Design notes that drive the code
--------------------------------
* Femoral neck DXA exists only in 2005-2006, 2007-2008, 2009-2010, 2013-2014 and
  2017-2018 (DXXFEM_*). 2011-2012 and 2015-2016 carry whole-body DXA only
  (DXX_*, no proximal femur), so they cannot contribute femoral neck BMD.
* Osteoporosis self-report (OSQ060) was not collected in 2011-2016.
* T-scores are computed against the sex-specific young-adult (20-29 y) femoral
  neck reference measured *within these same cycles* (D-F), so the reference is
  reproducible from the raw files instead of hard-coded literature constants.
* Linked mortality: fixed-width public-use 2019 LMF, positions confirmed from the
  CDC data dictionary and empirically (see module README). UCOD_LEADING 001 =
  diseases of heart, 005 = cerebrovascular disease; for 2015-2018 only 001/002/010
  are released, so cerebrovascular deaths are not separable in cycles I/J.

Outputs
  data/processed/nhanes/analysis_dataset.csv
  data/processed/nhanes/variable_dictionary.csv
  data/processed/nhanes/t_score_reference.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAW = Path("data/raw/nhanes")
OUT = Path("data/processed/nhanes")
OUT.mkdir(parents=True, exist_ok=True)

SUFFIX = {"2005-2006": "D", "2007-2008": "E", "2009-2010": "F", "2011-2012": "G",
          "2013-2014": "H", "2015-2016": "I", "2017-2018": "J"}
# cycles with femoral neck DXA
BMD_CYCLES = ["2005-2006", "2007-2008", "2009-2010", "2013-2014", "2017-2018"]
MORT_CYCLES = list(SUFFIX)

# fixed-width layout of NHANES_<cycle>_MORT_2019_PUBLIC.dat (1-based, inclusive)
LMF_COLS = {
    "SEQN": (1, 6), "ELIGSTAT": (15, 15), "MORTSTAT": (16, 16), "UCOD_LEADING": (17, 19),
    "DIABETES_MCD": (20, 20), "HYPERTEN_MCD": (21, 21), "DODQTR": (22, 22),
    "DODYEAR": (23, 26), "PERMTH_INT": (43, 45), "PERMTH_EXM": (46, 48),
}


def read(cycle, comp, columns=None):
    p = RAW / cycle / f"{comp}_{SUFFIX[cycle]}.XPT"
    if not p.exists():
        return None
    df = pd.read_sas(p, format="xport")
    if columns is not None:
        keep = [c for c in columns if c in df.columns]
        df = df[keep]
    df["SEQN"] = df["SEQN"].astype("int64")
    return df


def read_mortality(cycle):
    p = RAW / "mortality" / f"NHANES_{cycle.replace('-', '_')}_MORT_2019_PUBLIC.dat"
    if not p.exists():
        return None
    lines = p.read_text(errors="replace").splitlines()

    def cut(line, span):
        a, b = span
        v = line.ljust(48)[a - 1:b].strip()
        return v

    recs = []
    for ln in lines:
        if not ln.strip():
            continue
        r = {k: cut(ln, span) for k, span in LMF_COLS.items()}
        recs.append(r)
    df = pd.DataFrame(recs)
    for c in df.columns:
        df[c] = pd.to_numeric(df[c].replace({"": np.nan, ".": np.nan}), errors="coerce")
    df["SEQN"] = df["SEQN"].astype("int64")
    df["cycle"] = cycle
    return df


def egfr_ckdepi2021(scr, age, female):
    """CKD-EPI 2021 creatinine equation (race-free)."""
    k = np.where(female, 0.7, 0.9)
    a = np.where(female, -0.241, -0.302)
    ratio = scr / k
    return (142 * np.minimum(ratio, 1) ** a * np.maximum(ratio, 1) ** -1.200
            * 0.9938 ** age * np.where(female, 1.012, 1.0))


def build_cycle(cycle):
    demo = read(cycle, "DEMO", ["SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "DMDEDUC2",
                                "INDFMPIR", "WTMEC2YR", "SDMVPSU", "SDMVSTRA"])
    fem = read(cycle, "DXXFEM", ["SEQN", "DXXNKBMD", "DXXOFBMD", "DXAFMRST"])
    spn = read(cycle, "DXXSPN", ["SEQN", "DXXOSBMD"])
    parts = [demo]
    for d in (fem, spn, read(cycle, "BMX", ["SEQN", "BMXBMI", "BMXWAIST"]),
              read(cycle, "OSQ", ["SEQN", "OSQ060", "OSQ070"]),
              read(cycle, "MCQ", ["SEQN", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F"]),
              read(cycle, "BPQ", ["SEQN", "BPQ020", "BPQ040A", "BPQ050A"]),
              read(cycle, "BPX", ["SEQN", "BPXSY1", "BPXDI1"]),
              read(cycle, "DIQ", ["SEQN", "DIQ010"]),
              read(cycle, "SMQ", ["SEQN", "SMQ020", "SMQ040"]),
              read(cycle, "TCHOL", ["SEQN", "LBXTC"]),
              read(cycle, "HDL", ["SEQN", "LBDHDD"]),
              read(cycle, "BIOPRO", ["SEQN", "LBXSCR", "LBXSAL", "LBXSCA", "LBXSPH", "LBXSAPSI"]),
              read(cycle, "GHB", ["SEQN", "LBXGH"]),
              read(cycle, "RHQ", ["SEQN", "RHQ031", "RHQ060"])):
        if d is not None:
            d = d.drop_duplicates("SEQN")
            parts.append(d)
    df = parts[0]
    for p in parts[1:]:
        df = df.merge(p, on="SEQN", how="left")
    df["cycle"] = cycle
    return df


def main():
    frames = [build_cycle(c) for c in BMD_CYCLES]
    raw = pd.concat(frames, ignore_index=True)
    print(f"stacked BMD cycles: {raw.shape}")

    # ---- T-score reference: sex-specific young adults (20-29 y) in cycles D-F ----
    ya = raw[(raw["RIDAGEYR"].between(20, 29)) & raw["DXXNKBMD"].notna()
             & raw["cycle"].isin(["2005-2006", "2007-2008", "2009-2010"])]
    ref = (ya.groupby("RIAGENDR")["DXXNKBMD"]
             .agg(n="count", mean="mean", sd="std").reset_index())
    ref["sex"] = ref["RIAGENDR"].map({1: "Male", 2: "Female"})
    # combine men and women for a sex-neutral fallback reference
    allref = pd.DataFrame([{"n": len(ya), "mean": ya["DXXNKBMD"].mean(),
                            "sd": ya["DXXNKBMD"].std(), "sex": "All"}])
    ref = pd.concat([ref[["sex", "n", "mean", "sd"]], allref], ignore_index=True)
    ref.to_csv(OUT / "t_score_reference.csv", index=False, float_format="%.6f")
    print(ref.to_string(index=False))

    mean_map = dict(zip(ref["sex"], ref["mean"]))
    sd_map = dict(zip(ref["sex"], ref["sd"]))
    raw["sex"] = raw["RIAGENDR"].map({1: "Male", 2: "Female"})
    raw["t_score_fn"] = (raw["DXXNKBMD"] - raw["sex"].map(mean_map)) / raw["sex"].map(sd_map)
    raw["t_score_all"] = (raw["DXXNKBMD"] - mean_map["All"]) / sd_map["All"]

    # ---- derived variables ----
    d = pd.DataFrame({"SEQN": raw["SEQN"], "cycle": raw["cycle"], "sex": raw["sex"],
                      "age": raw["RIDAGEYR"], "race": raw["RIDRETH1"].map(
                          {1: "Mexican American", 2: "Other Hispanic", 3: "NH White",
                           4: "NH Black", 5: "Other/Multi"}),
                      "education": pd.cut(raw["DMDEDUC2"], [0, 2, 3, 5],
                                          labels=["<High school", "High school/GED", ">High school"]),
                      "pir": raw["INDFMPIR"], "bmi": raw["BMXBMI"], "waist": raw["BMXWAIST"],
                      "fn_bmd": raw["DXXNKBMD"], "tf_bmd": raw["DXXOFBMD"], "ls_bmd": raw["DXXOSBMD"],
                      "t_score": raw["t_score_fn"], "t_score_all": raw["t_score_all"],
                      "wt_mec": raw["WTMEC2YR"], "psu": raw["SDMVPSU"], "strata": raw["SDMVSTRA"]})

    # cardiovascular disease (self-reported physician diagnosis)
    cvd_items = {"mcq_chf": "MCQ160B", "mcq_chd": "MCQ160C", "mcq_angina": "MCQ160D",
                 "mcq_mi": "MCQ160E", "mcq_stroke": "MCQ160F"}
    for name, col in cvd_items.items():
        d[name] = (raw[col] == 1).where(raw[col].isin([1, 2]))
    answered = raw[list(cvd_items.values())].isin([1, 2]).any(axis=1)
    d["cvd"] = d[list(cvd_items)].any(axis=1).where(answered)
    d["chd"] = (d["mcq_chd"] | d["mcq_mi"]).where(d["mcq_chd"].notna() | d["mcq_mi"].notna())
    d["stroke"] = d["mcq_stroke"]
    d["hf"] = d["mcq_chf"]

    # osteoporosis
    d["self_op"] = (raw["OSQ060"] == 1).where(raw["OSQ060"].isin([1, 2]))
    d["menopause"] = (raw["RHQ031"] == 2).where(raw["RHQ031"].isin([1, 2]) & (d["sex"] == "Female"))
    d["age_menopause"] = raw["RHQ060"].where(raw["RHQ060"] < 60)
    # WHO categories on the sex-specific T-score
    d["bmd_cat"] = pd.cut(d["t_score"], [-np.inf, -2.5, -1.0, np.inf],
                          labels=["Osteoporosis", "Osteopenia", "Normal"])
    d["bmd_cat_all"] = pd.cut(d["t_score_all"], [-np.inf, -2.5, -1.0, np.inf],
                              labels=["Osteoporosis", "Osteopenia", "Normal"])

    # hypertension: self-report OR measured SBP>=140 / DBP>=90 OR medication
    measured_htn = (raw["BPXSY1"] >= 140) | (raw["BPXDI1"] >= 90)
    d["hypertension"] = (raw["BPQ020"].eq(1) | measured_htn | raw["BPQ050A"].eq(1)) \
        .where(raw["BPQ020"].isin([1, 2]) | measured_htn.notna())
    # diabetes: self-report OR HbA1c >= 6.5 %
    d["diabetes"] = (raw["DIQ010"].eq(1) | (raw["LBXGH"] >= 6.5)) \
        .where(raw["DIQ010"].isin([1, 2]) | raw["LBXGH"].notna())
    # smoking
    d["smoking"] = np.select(
        [raw["SMQ020"].eq(2), raw["SMQ020"].eq(1) & raw["SMQ040"].eq(3),
         raw["SMQ020"].eq(1) & raw["SMQ040"].isin([1, 2])],
        ["Never", "Former", "Current"], default=None)
    d["total_chol"] = raw["LBXTC"]
    d["hdl"] = raw["LBDHDD"]
    d["creatinine"] = raw["LBXSCR"]
    d["albumin"] = raw["LBXSAL"]
    d["calcium"] = raw["LBXSCA"]
    d["phosphate"] = raw["LBXSPH"]
    d["alp"] = raw["LBXSAPSI"]
    d["hba1c"] = raw["LBXGH"]
    d["egfr"] = egfr_ckdepi2021(raw["LBXSCR"].to_numpy(float), raw["RIDAGEYR"].to_numpy(float),
                                (raw["RIAGENDR"] == 2).to_numpy())
    d.loc[raw["LBXSCR"].isna(), "egfr"] = np.nan

    # ---- mortality ----
    mort = pd.concat([read_mortality(c) for c in MORT_CYCLES], ignore_index=True)
    mort = mort.drop_duplicates(["SEQN", "cycle"])
    mort["heart_death"] = (mort["UCOD_LEADING"] == 1).astype(float).where(mort["MORTSTAT"] == 1, 0.0)
    mort["cerebro_death"] = (mort["UCOD_LEADING"] == 5).astype(float).where(mort["MORTSTAT"] == 1, 0.0)
    mort["cvd_death"] = (mort["MORTSTAT"] == 1) & mort["UCOD_LEADING"].isin([1, 5])
    d = d.merge(mort[["SEQN", "cycle", "ELIGSTAT", "MORTSTAT", "UCOD_LEADING",
                      "cvd_death", "heart_death", "cerebro_death",
                      "PERMTH_INT", "PERMTH_EXM"]], on=["SEQN", "cycle"], how="left")
    d["fu_years_exm"] = d["PERMTH_EXM"] / 12.0

    # ---- encode logical columns as 0/1 integers so downstream R keeps them numeric
    bool_cols = ["cvd", "chd", "stroke", "hf", "self_op", "menopause", "hypertension",
                 "diabetes", "cvd_death", "heart_death", "cerebro_death"]
    for c in bool_cols:
        d[c] = d[c].map({True: 1, False: 0}).astype("Int64")

    # ---- analytic sample ----
    n0 = len(d)
    d = d[d["age"] >= 40]
    n1 = len(d)
    d = d[d["fn_bmd"].notna() & d["cvd"].notna() & (d["ELIGSTAT"] == 1)]
    print(f"rows: {n0} -> age>=40: {n1} -> with femur BMD + CVD + mortality-eligible: {len(d)}")
    print(d.groupby("cycle").size().to_string())

    d.to_csv(OUT / "analysis_dataset.csv", index=False, float_format="%.6f")

    vardict = [
        ("SEQN", "NHANES respondent sequence number"), ("cycle", "2-year survey cycle"),
        ("sex", "RIAGENDR"), ("age", "RIDAGEYR (years)"), ("race", "RIDRETH1"),
        ("education", "DMDEDUC2 collapsed"), ("pir", "INDFMPIR poverty income ratio"),
        ("bmi", "BMXBMI"), ("waist", "BMXWAIST"),
        ("fn_bmd", "DXXNKBMD femoral neck BMD (g/cm2), primary exposure"),
        ("tf_bmd", "DXXOFBMD total femur BMD"), ("ls_bmd", "DXXOSBMD total lumbar spine BMD"),
        ("t_score", "sex-specific T-score vs young-adult reference of this dataset"),
        ("t_score_all", "sex-neutral T-score (sensitivity)"),
        ("bmd_cat", "WHO category from t_score"), ("self_op", "OSQ060 self-reported osteoporosis"),
        ("cvd", "any of CHF/CHD/angina/MI/stroke, self-reported physician diagnosis"),
        ("chd", "coronary heart disease or MI"), ("stroke", "self-reported stroke"),
        ("hf", "self-reported congestive heart failure"),
        ("hypertension", "self-report or SBP>=140 or DBP>=90 or medication"),
        ("diabetes", "self-report or HbA1c>=6.5%"), ("smoking", "never/former/current"),
        ("total_chol", "LBXTC mg/dL"), ("hdl", "LBDHDD mg/dL"),
        ("albumin", "LBXSAL g/dL"), ("calcium", "LBXSCA mg/dL"), ("phosphate", "LBXSPH mg/dL"),
        ("alp", "LBXSAPSI alkaline phosphatase IU/L"), ("hba1c", "LBXGH %"),
        ("egfr", "CKD-EPI 2021 creatinine eGFR mL/min/1.73m2"),
        ("menopause", "no period in past 12 months (women)"),
        ("MORTSTAT", "0 assumed alive, 1 assumed deceased (through 2019-12-31)"),
        ("UCOD_LEADING", "1 heart disease, 2 cancer, 3 CLRD, 4 accidents, 5 cerebrovascular"),
        ("cvd_death", "death from heart disease (001) or cerebrovascular disease (005)"),
        ("fu_years_exm", "follow-up years from MEC exam (PERMTH_EXM/12)"),
        ("wt_mec", "WTMEC2YR 2-year MEC weight; divide by 5 when pooling 5 cycles"),
        ("psu", "SDMVPSU"), ("strata", "SDMVSTRA"),
    ]
    pd.DataFrame(vardict, columns=["variable", "description"]).to_csv(
        OUT / "variable_dictionary.csv", index=False)
    print(f"-> {OUT/'analysis_dataset.csv'}")


if __name__ == "__main__":
    sys.exit(main())

