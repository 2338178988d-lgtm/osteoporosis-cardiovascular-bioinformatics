"""Shared reader for GWAS Catalog harmonised summary statistics.

The six files this project uses do not all carry the same column names:

  eBMD (Loh 2018)        hm_rsid, hm_chrom, hm_pos, hm_beta, hm_effect_allele, ..., n
  MI/STROKE/HF (Sakaue 2021, Mishra 2022, Levin 2022)
                         rsid, chromosome, base_pair_location, beta, standard_error,
                         effect_allele, other_allele, effect_allele_frequency, p_value

so every file is mapped onto one schema (SNP, chr, pos, effect_allele, other_allele,
eaf, beta, se, p, n) on read. Where a file has no sample-size column the analysis N is
taken from `TRAIT_N`, which records the value and where it came from; these are used only
for the Steiger test and to compute instrument R^2, never to derive an effect estimate.
"""
import gzip
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data" / "raw" / "gwas"
FINNGEN = ROOT / "data" / "raw" / "finngen"

# a key whose file lives outside data/raw/gwas/ is listed here with its own directory
DIRS = {}

FILES = {
    "eBMD": "GCST90029004_eBMD.h.tsv.gz",
    "OP": "GCST90018887_OP.h.tsv.gz",
    "CAD": "GCST90132315_CAD.h.tsv.gz",
    "MI": "GCST90018877_MI.h.tsv.gz",
    "STROKE": "GCST90104534_STROKE.h.tsv.gz",
    "HF": "GCST90162626_HF.h.tsv.gz",
    # added for the eBMD -> heart-failure follow-up (module 05)
    "HF_SHAH": "HF_SHAH_GCST009541.h.tsv.gz",
    "HF_CHRONIC": "HF_CHRONIC_GCST90018806.h.tsv.gz",
    "BMI": "BMI_GCST90029007.h.tsv.gz",
    "SBP": "SBP_GCST90000066.h.tsv.gz",
    "T2D": "T2D_GCST005413.h.tsv.gz",
    # bone-derived mediators (module 06)
    "FGF23": "FGF23_GCST90469223.h.tsv.gz",
    "SOST": "SOST_GCST90470711.h.tsv.gz",
    "BGLAP": "BGLAP_GCST90468440.h.tsv.gz",
    "CST3": "CST3_GCST90025945.h.tsv.gz",
    "CRP": "CRP_GCST90025959.h.tsv.gz",
    # extended mediators and independent populations (module 08)
    "DKK1": "DKK1_GCST90468997.h.tsv.gz",
    "PTH": "PTH_GCST90470372.h.tsv.gz",
    "TESTO": "TESTO_GCST90012114.h.tsv.gz",
    "ESTRADIOL": "ESTRADIOL_GCST90012105.h.tsv.gz",
    "VITD": "VITD_GCST90000618.h.tsv.gz",
    "FG_EA_HF": "FG_EA_HF_GCST90018586.h.tsv.gz",
    "FG_EA_OP": "FG_EA_OP_GCST90018667.h.tsv.gz",
    "HFPEUR": "HFPEUR_GCST90726619.h.tsv.gz",
    "HFREUR": "HFREUR_GCST90726620.h.tsv.gz",
    "HFPEA": "HFPEA_GCST90668011.h.tsv.gz",
    "HFREA": "HFREA_GCST90668010.h.tsv.gz",
    # FinnGen R12: independent population, no UK Biobank overlap
    "FG_HF": "finngen_R12_I9_HEARTFAIL.gz",
    "FG_OP": "finngen_R12_M13_OSTEOPOROSIS.gz",
    "FG_HF_OBESE": "finngen_R12_I9_HEARTFAIL_AND_OVERWEIGHT.gz",
    "FG_HF_CHD": "finngen_R12_I9_HEARTFAIL_AND_CHD.gz",
    # body-composition alternative explanation for the eBMD effect (module 11)
    "ALM": "ALM_GCST90000025.h.tsv.gz",
    "HEIGHT": "HEIGHT_GCST90029008.h.tsv.gz",
    # exposure-robustness work (module 13)
    "FNBMD": "FNBMD_GCST90446641.h.tsv.gz",
    "LVM": "LVM_GCST90244711.h.tsv.gz",
    "LVMI": "LVMI_GCST90435260.h.tsv.gz",
    # FinnGen negative-control endpoints (module 13)
    "FN_NEG_APPENDACUT": "finngen_R12_K11_APPENDACUT.gz",
    "FN_NEG_CATARACT": "finngen_R12_H7_CATARACTSENILE.gz",
    "FN_NEG_ATOPIC": "finngen_R12_L12_ATOPIC.gz",
    "FN_NEG_MIGRAINE": "finngen_R12_G6_MIGRAINE.gz",
    "FN_NEG_TONSILL": "finngen_R12_J10_TONSILLECTOMY.gz",
    "FN_NEG_SINUSITIS": "finngen_R12_J10_SINUSITIS.gz",
    # cardiac ladder (module 14): FinnGen R12, read through the remote tabix reader,
    # never materialised as a full file
    "FG_AF": "finngen_R12_I9_AF.gz",
    "FG_HYPERTENS": "finngen_R12_I9_HYPTENS.gz",
    "FG_HYPERTENS_HD": "finngen_R12_I9_HYPTENSHD.gz",
    "FG_NONISCH_CM": "finngen_R12_I9_NONISCHCARDMYOP.gz",
    "FG_CARDMYO": "finngen_R12_I9_CARDMYO.gz",
    "FG_CAVS": "finngen_R12_I9_CAVS_OPERATED.gz",
    "FG_VHD": "finngen_R12_I9_VHD.gz",
    "FG_MI_STRICT": "finngen_R12_I9_MI_STRICT.gz",
    # discriminating panel + dilution series (module 14 extension)
    "FG_VTE": "finngen_R12_I9_VTE.gz",
    "FG_VARICVE": "finngen_R12_I9_VARICVE.gz",
    "FG_AORTANEUR": "finngen_R12_I9_AORTANEUR.gz",
    "FG_PULMEMB": "finngen_R12_I9_PULMEMB.gz",
    "FG_OTHARR": "finngen_R12_I9_OTHARR.gz",
    "FG_AF_REIMB": "finngen_R12_I9_AF_REIMB.gz",
    "FG_CONDUCTIO": "finngen_R12_I9_CONDUCTIO.gz",
    "FG_AVBLOCK": "finngen_R12_I9_AVBLOCK.gz",
    "FG_PAROXTAC": "finngen_R12_I9_PAROXTAC.gz",
    "FG_CVD_ALL": "finngen_R12_FG_CVD.gz",
    "FG_IHD": "finngen_R12_I9_IHD.gz",
    "FG_CHD": "finngen_R12_I9_CHD.gz",
    # aortic valve calcification (module 14 extension, non-Finnish control)
    "AVC": "AVC_GCST90310293.h.tsv.gz",
    # modality triangulation (module 15)
    "FNBMD_2": "FNBMD_2_GCST90446628.h.tsv.gz",
    "FNBMD_3": "FNBMD_3_GCST90568458.h.tsv.gz",
    "QUS_2": "QUS_2_GCST90025982.h.tsv.gz",
    "QUS_3": "QUS_3_GCST006979.h.tsv.gz",
}
for _k in ("FG_HF", "FG_OP", "FG_HF_OBESE", "FG_HF_CHD", "FN_NEG_APPENDACUT",
           "FN_NEG_CATARACT", "FN_NEG_ATOPIC", "FN_NEG_MIGRAINE", "FN_NEG_TONSILL",
           "FN_NEG_SINUSITIS", "FG_AF", "FG_HYPERTENS", "FG_HYPERTENS_HD",
           "FG_NONISCH_CM", "FG_CARDMYO", "FG_CAVS", "FG_VHD", "FG_MI_STRICT",
           "FG_VTE", "FG_VARICVE", "FG_AORTANEUR", "FG_PULMEMB", "FG_OTHARR",
           "FG_AF_REIMB", "FG_CONDUCTIO", "FG_AVBLOCK", "FG_PAROXTAC", "FG_CVD_ALL",
           "FG_IHD", "FG_CHD"):
    DIRS[_k] = FINNGEN

# fallback analysis sample size, and where it came from
TRAIT_N = {
    "eBMD": (445921, "n column in the file"),
    "OP": (479981, "GCST90018887 study record, Sakaue 2021"),
    "CAD": (1167328, "GCST90132315 study record, Aragam 2022"),
    "MI": (440906, "GCST90018877 study record, Sakaue 2021"),
    "STROKE": (1234808, "GCST90104534 study record, Mishra 2022"),
    "HF": (1270968, "GCST90162626 study record, Levin 2022"),
    "HF_SHAH": (977323, "GCST009541, Shah 2020 (47,309 cases / 930,014 controls)"),
    "HF_CHRONIC": (486160, "GCST90018806, Sakaue 2021 (14,262 cases / 471,898 controls)"),
    "BMI": (532396, "GCST90029007, Loh 2018"),
    "SBP": (810865, "GCST90000066, Surendran 2020"),
    "T2D": (898130, "GCST005413, Mahajan 2018 (74,124 cases / 824,006 controls)"),
    "FGF23": (47745, "GCST90469223, Loya 2025"),
    "SOST": (47745, "GCST90470711, Loya 2025"),
    "BGLAP": (47745, "GCST90468440, Loya 2025"),
    "CST3": (437846, "GCST90025945, Barton 2021"),
    "CRP": (436939, "GCST90025959, Barton 2021"),
    "DKK1": (47745, "GCST90468997, Loya 2025"),
    "PTH": (47745, "GCST90470372, Loya 2025"),
    "TESTO": (425097, "GCST90012114, Ruth 2020"),
    "ESTRADIOL": (206927, "GCST90012105, Ruth 2020"),
    "VITD": (496946, "GCST90000618, vitamin D"),
    "FG_EA_HF": (178726, "GCST90018586, Sakaue 2021 BioBank Japan"),
    "FG_EA_OP": (178726, "GCST90018667, Sakaue 2021 BioBank Japan"),
    "HFPEUR": (756189, "GCST90726619, Henry 2025 (3,590 HFpEF / 752,599)"),
    "HFREUR": (940805, "GCST90726620, Henry 2025 (4,975 HFrEF / 935,830)"),
    "HFPEA": (201831, "GCST90668011, Enzan 2025 East Asian HFpEF"),
    "HFREA": (201831, "GCST90668010, Enzan 2025 East Asian HFrEF"),
    "FG_HF": (500348, "FinnGen R12 manifest: 37,653 cases / 462,695 controls"),
    "FG_OP": (483725, "FinnGen R12 manifest: 10,461 cases / 473,264 controls"),
    "FG_HF_OBESE": (229541, "FinnGen R12 manifest: 22,965 cases / 206,576 controls"),
    "FG_HF_CHD": (461230, "FinnGen R12 manifest: 17,532 cases / 443,698 controls"),
    "ALM": (450243, "GCST90000025, Pei 2020 (appendicular lean mass, UKB)"),
    "HEIGHT": (673878, "GCST90029008, Loh 2018 (standing height, UKB)"),
    "FNBMD": (79940, "GCST90446641 (DXA femoral neck BMD, European)"),
    "LVM": (42432, "GCST90244711 (left ventricular mass, 41,855 EUR + 577 Asian)"),
    "LVMI": (36083, "GCST90435260 (LV mass indexed to body surface area)"),
    "FN_NEG_APPENDACUT": (497357, "FinnGen R12 manifest: 38,933 cases / 458,424 controls"),
    "FN_NEG_CATARACT": (493421, "FinnGen R12 manifest: 83,886 cases / 409,535 controls"),
    "FN_NEG_ATOPIC": (464119, "FinnGen R12 manifest: 31,245 cases / 432,874 controls"),
    "FN_NEG_MIGRAINE": (401499, "FinnGen R12 manifest: 26,894 cases / 374,605 controls"),
    "FN_NEG_TONSILL": (500348, "FinnGen R12 manifest: 39,135 cases / 461,213 controls"),
    "FN_NEG_SINUSITIS": (431268, "FinnGen R12 manifest: 28,616 cases / 402,652 controls"),
    # module 14, case/control counts from the FinnGen R12 manifest
    "FG_AF": (316342, "FinnGen R12 manifest: 63,532 cases / 252,810 controls"),
    "FG_HYPERTENS": (500264, "FinnGen R12 manifest: 154,630 cases / 345,634 controls"),
    "FG_HYPERTENS_HD": (357362, "FinnGen R12 manifest: 11,728 cases / 345,634 controls"),
    "FG_NONISCH_CM": (416260, "FinnGen R12 manifest: 13,917 cases / 402,343 controls"),
    "FG_CARDMYO": (382992, "FinnGen R12 manifest: 7,649 cases / 375,343 controls"),
    "FG_CAVS": (500348, "FinnGen R12 manifest: 12,418 cases / 487,930 controls"),
    "FG_VHD": (403945, "FinnGen R12 manifest: 28,602 cases / 375,343 controls"),
    "FG_MI_STRICT": (447837, "FinnGen R12 manifest: 31,666 cases / 416,171 controls"),
    # module 14 extension, case/control counts from the FinnGen R12 manifest
    "FG_VTE": (500348, "FinnGen R12 manifest: 26,333 cases / 474,015 controls"),
    "FG_VARICVE": (470690, "FinnGen R12 manifest: 38,467 cases / 432,223 controls"),
    "FG_AORTANEUR": (473043, "FinnGen R12 manifest: 9,937 cases / 463,106 controls"),
    "FG_PULMEMB": (499081, "FinnGen R12 manifest: 12,762 cases / 486,319 controls"),
    "FG_OTHARR": (414608, "FinnGen R12 manifest: 39,265 cases / 375,343 controls"),
    "FG_AF_REIMB": (293061, "FinnGen R12 manifest: 40,251 cases / 252,810 controls"),
    "FG_CONDUCTIO": (389370, "FinnGen R12 manifest: 14,027 cases / 375,343 controls"),
    "FG_AVBLOCK": (383193, "FinnGen R12 manifest: 7,850 cases / 375,343 controls"),
    "FG_PAROXTAC": (265688, "FinnGen R12 manifest: 12,878 cases / 252,810 controls"),
    "FG_CVD_ALL": (500348, "FinnGen R12 manifest: 247,538 cases / 252,810 controls"),
    "FG_IHD": (500348, "FinnGen R12 manifest: 84,088 cases / 416,260 controls"),
    "FG_CHD": (500348, "FinnGen R12 manifest: 56,650 cases / 443,698 controls"),
    "AVC": (956682, "GCST90310293 (aortic valve calcification, 14,819 / 941,863)"),
    # module 15
    "FNBMD_2": (47205, "GCST90446628 (DXA femoral neck BMD, European)"),
    "FNBMD_3": (32017, "GCST90568458 (DXA femoral neck BMD, European)"),
    "QUS_2": (445855, "GCST90025982 (heel ultrasound BMD T-score, European)"),
    "QUS_3": (426824, "GCST006979 (heel ultrasound BMD, British)"),
}

# case / control split for case-control traits, used only for power statements
CASE_COUNTS = {
    "HF": (47309, 930014),
    "FG_HF": (37653, 462695),
    "FG_HF_OBESE": (22965, 206576),
    "FG_HF_CHD": (17532, 443698),
    "FN_NEG_APPENDACUT": (38933, 458424),
    "FN_NEG_CATARACT": (83886, 409535),
    "FN_NEG_ATOPIC": (31245, 432874),
    "FN_NEG_MIGRAINE": (26894, 374605),
    "FN_NEG_TONSILL": (39135, 461213),
    "FN_NEG_SINUSITIS": (28616, 402652),
    # cardiac ladder (module 14)
    "FG_AF": (63532, 252810),
    "FG_HYPERTENS": (154630, 345634),
    "FG_HYPERTENS_HD": (11728, 345634),
    "FG_NONISCH_CM": (13917, 402343),
    "FG_CARDMYO": (7649, 375343),
    "FG_CAVS": (12418, 487930),
    "FG_VHD": (28602, 375343),
    "FG_MI_STRICT": (31666, 416171),
    # module 14 extension
    "FG_VTE": (26333, 474015),
    "FG_VARICVE": (38467, 432223),
    "FG_AORTANEUR": (9937, 463106),
    "FG_PULMEMB": (12762, 486319),
    "FG_OTHARR": (39265, 375343),
    "FG_AF_REIMB": (40251, 252810),
    "FG_CONDUCTIO": (14027, 375343),
    "FG_AVBLOCK": (7850, 375343),
    "FG_PAROXTAC": (12878, 252810),
    "FG_CVD_ALL": (247538, 252810),
    "FG_IHD": (84088, 416260),
    "FG_CHD": (56650, 443698),
    "AVC": (14819, 941863),
}

# (target name, candidate source columns in order of preference)
SCHEMA = [
    ("SNP", ["hm_rsid", "rsid", "variant_id", "rsids"]),    ("chr", ["hm_chrom", "chromosome", "#chrom"]),
    ("pos", ["hm_pos", "base_pair_location", "pos"]),
    ("effect_allele", ["hm_effect_allele", "effect_allele", "alt"]),
    ("other_allele", ["hm_other_allele", "other_allele", "ref"]),
    ("eaf", ["hm_effect_allele_frequency", "effect_allele_frequency", "af_alt"]),
    ("beta", ["hm_beta", "beta"]),
    ("se", ["standard_error", "sebeta"]),
    ("p", ["p_value", "pval"]),
    ("n", ["n", "N"]),
]


# Traits whose file carries an `n` column that is not a sample size. GCST006979 ships a
# column named `n` whose values sit between 0 and 1 (an allele-frequency-like quantity); using
# it as N makes the instrument R^2 collapse to ~1 and reports an impossible Sum R2 of 526.
# These keys fall back to TRAIT_N instead.
BAD_N = {"QUS_3"}


def path_of(key):
    return DIRS.get(key, RAW) / FILES[key]


def header(key):
    with gzip.open(path_of(key), "rt") as f:
        return f.readline().rstrip("\n").split("\t")


def column_map(key):
    h = set(header(key))
    if key in BAD_N:
        h.discard("n")
        h.discard("N")
    mapping, missing = {}, []
    for target, cands in SCHEMA:
        src = next((c for c in cands if c in h), None)
        if src is None:
            missing.append(target)
        else:
            mapping[src] = target
    if "SNP" in missing or "beta" in missing or "se" in missing or "p" in missing:
        raise RuntimeError(f"{key}: required columns missing from {FILES[key]} ({missing})")
    return mapping, missing


def read_rows(key, p_thresh, chunksize=2_000_000):
    """Stream the file, returning only rows with p < p_thresh on the common schema."""
    mapping, missing = column_map(key)
    n_fallback = TRAIT_N[key][0]
    keep = []
    for chunk in pd.read_csv(path_of(key), sep="\t", usecols=list(mapping), chunksize=chunksize,
                             low_memory=False, dtype={"rsid": str, "hm_rsid": str, "variant_id": str, "rsids": str}):
        d = chunk.rename(columns=mapping)
        d = d[pd.to_numeric(d["p"], errors="coerce") < p_thresh]
        if len(d):
            keep.append(d)
    if not keep:
        empty = pd.DataFrame({t: pd.Series(dtype="float64") for t, _ in SCHEMA})
        empty["SNP"] = pd.Series(dtype="object")
        return empty
    d = pd.concat(keep, ignore_index=True)
    if "n" in d:
        d["n"] = pd.to_numeric(d["n"], errors="coerce").fillna(n_fallback)
    else:
        d["n"] = float(n_fallback)
    return d


def read_at_snps(key, snps, chunksize=2_000_000):
    """Stream the file once, keeping the requested SNPs on the common schema."""
    mapping, _ = column_map(key)
    n_fallback = TRAIT_N[key][0]
    want = set(snps)
    rows = []
    for chunk in pd.read_csv(path_of(key), sep="\t", usecols=list(mapping), chunksize=chunksize,
                             low_memory=False, dtype={"rsid": str, "hm_rsid": str, "variant_id": str, "rsids": str}):
        chunk = chunk.rename(columns=mapping)
        hit = chunk[chunk["SNP"].isin(want)]
        if len(hit):
            rows.append(hit)
    if not rows:
        return None
    d = pd.concat(rows, ignore_index=True)
    if "n" in d:
        d["n"] = pd.to_numeric(d["n"], errors="coerce").fillna(n_fallback)
    else:
        d["n"] = float(n_fallback)
    return d.dropna(subset=["beta", "se"]).drop_duplicates("SNP")

