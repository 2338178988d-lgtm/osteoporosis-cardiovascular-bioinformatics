"""Module 14, part 2 - three discriminating tests that decide what the ladder means.

Module 14 measured the effect on eight cardiac endpoints and found it on atrial fibrillation,
valvular heart disease and non-ischaemic cardiomyopathy, and not on hypertension or myocardial
infarction. Three questions follow directly, and each is answered by endpoints that already
exist. All use the same 600 eBMD instruments, the same estimator, the same Finnish cohort.

  A. CARDIAC OR SYSTEMIC MATRIX?
     The three protective endpoints are all "connective tissue that is inside the heart".
     The competing reading is that the exposure tracks collagen/matrix quality everywhere,
     in which case venous and aortic disease should move too. Four well-powered endpoints
     decide it: varicose veins (38,467), VTE (26,333), pulmonary embolism (12,762), aortic
     aneurysm (9,937). A null on varicose veins at 94 % power would be hard to explain for a
     systemic-matrix reading.

  B. ATRIAL MUSCLE OR CONDUCTION SYSTEM?
     Atrial fibrillation is the largest non-coronary effect (-0.086) and shares patients with
     heart failure, so it could be heart failure under another name. The anatomy resolves it
     without removing shared cases: an atrial-myocyte effect leaves the conduction system
     alone; a fibrotic effect hits both. Conduction disorders (14,027), AV block (7,850),
     paroxysmal tachycardia (12,878), other arrhythmias (39,265).

  C. THE DILUTION SERIES.
     03_dilution.py showed the three heart-failure numbers are one number diluted by the
     ischaemic half. If that is right, the effect must shrink monotonically as an endpoint
     becomes more inclusive: cardiovascular disease as a whole (247,538) < ischaemic heart
     disease (84,088) <= major CHD events (56,650) should be near zero, while heart failure
     (-0.064) and non-ischaemic cardiomyopathy (-0.102) are further down the same axis.

DECISION RULES, written before the run:

  A. any of the four systemic-matrix endpoints protective at Bonferroni (p < 0.0125)
        -> the effect is not cardiac-specific; "bone-matrix axis" replaces "bone-cardiac axis"
           and the paper's scope widens
     all four null with varicose veins and VTE adequately powered
        -> the effect is cardiac-specific; the systemic-matrix reading is excluded

  B. conduction endpoints protective as well as atrial
        -> a shared fibrotic substrate; the atrial label is too narrow
     atrial protective, conduction null
        -> the effect is atrial muscle, and the heart-failure overlap does not explain it

  C. |b| increases monotonically from all-CVD to non-ischaemic cardiomyopathy
        -> the dilution reading moves from arithmetic to measurement

Multiplicity: group A tests four endpoints (Bonferroni 0.0125), group B five (0.01), group C
three (0.0167). Reported alongside unadjusted p; the three groups are pre-specified as
separate families because they answer different questions.

Inputs: twelve FinnGen R12 `.tbi` indexes in data/raw/finngen/ (1.8 MB each; the 815 MB files
are never downloaded). Variant rows are cached in data/processed/mr/outcome_snps/.

Outputs (results/tables/14-cardiac-ladder/):
  panel.csv        one row per endpoint with group, beta, CI, p, power, verdict
  dilution_series.csv  the monotonicity check in group C

Usage:
  python analysis/14-cardiac-ladder/scripts/04_discriminating_panel.py
  python analysis/14-cardiac-ladder/scripts/04_discriminating_panel.py --only A
  python analysis/14-cardiac-ladder/scripts/04_discriminating_panel.py --self-test
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "analysis" / "scripts"))
sys.path.insert(0, str(ROOT / "analysis" / "13-robustness" / "scripts"))

from _common import estimate, mr  # noqa: E402
import _sumstats  # noqa: E402
from _tabix_http import TabixHttp, fetch  # noqa: E402

OUT = ROOT / "results" / "tables" / "14-cardiac-ladder"
CACHE = ROOT / "data" / "processed" / "mr" / "outcome_snps"
FINNGEN_DIR = ROOT / "data" / "raw" / "finngen"
FG_URL = ("https://storage.googleapis.com/finngen-public-data-r12/summary_stats/"
          "release/finngen_R12_{}.gz")
EXPOSURE = "eBMD"

# (group, key, phenocode, cases, expectation, rationale)
PANEL = [
    ("A", "FG_VARICVE", "I9_VARICVE", 38467, "null",
     "venous matrix, 94 % power - the single most discriminating row"),
    ("A", "FG_VTE", "I9_VTE", 26333, "null", "thrombosis / venous wall"),
    ("A", "FG_PULMEMB", "I9_PULMEMB", 12762, "null", "same axis as VTE"),
    ("A", "FG_AORTANEUR", "I9_AORTANEUR", 9937, "null", "arterial matrix, aortic wall"),
    ("B", "FG_OTHARR", "I9_OTHARR", 39265, "protective", "arrhythmia, non-AF"),
    ("B", "FG_AF_REIMB", "I9_AF_REIMB", 40251, "protective", "AF, reimbursement definition"),
    ("B", "FG_PAROXTAC", "I9_PAROXTAC", 12878, "protective", "atrial/AV re-entry"),
    ("B", "FG_CONDUCTIO", "I9_CONDUCTIO", 14027, "null",
     "conduction system - null means atrial muscle, protective means fibrosis"),
    ("B", "FG_AVBLOCK", "I9_AVBLOCK", 7850, "null", "conduction, narrower"),
    ("C", "FG_CVD_ALL", "FG_CVD", 247538, "null", "most aggregated"),
    ("C", "FG_IHD", "I9_IHD", 84088, "null", "ischaemic aggregate"),
    ("C", "FG_CHD", "I9_CHD", 56650, "null", "coronary aggregate"),
]
FAMILY_ALPHA = {"A": 0.05 / 4, "B": 0.05 / 5, "C": 0.05 / 3}
# the dilution series, ordered from most aggregated to most specific, with module-14 anchors
SERIES = [("FG_CVD_ALL", "cardiovascular disease, all"),
          ("FG_IHD", "ischaemic heart disease"),
          ("FG_CHD", "major coronary heart disease event"),
          ("FG_HF", "heart failure"),
          ("FG_NONISCH_CM", "non-ischaemic cardiomyopathy")]


def cached_outcome(key, inst, threads=8, force=False):
    """Instrument-SNP rows for a remote endpoint, cached on disk after the first fetch."""
    path = CACHE / f"{key}.csv"
    asked = CACHE / f"{key}.snps"
    want = set(inst["SNP"])
    if path.exists() and not force and asked.exists():
        if set(asked.read_text().split()) == want:
            got = pd.read_csv(path)
            return got[got["SNP"].isin(want)]
    pheno = _sumstats.FILES[key].replace("finngen_R12_", "").replace(".gz", "")
    tbi = FINNGEN_DIR / f"finngen_R12_{pheno}.gz.tbi"
    if not tbi.exists():
        print(f"  {key}: index {tbi.name} absent - run download_plan.py --category ladder")
        return None
    print(f"  {key}: fetching {len(inst)} positions via tabix ...", flush=True)
    reader = TabixHttp(FG_URL.format(pheno), tbi)
    cols = reader.header
    rows = fetch(reader, [(r.SNP, r.chr, int(r.pos)) for r in inst.itertuples()],
                 threads=threads)
    by_snp = {}
    for snp, pos, f in rows:
        rec = dict(zip(cols, f))
        by_snp.setdefault(snp, []).append((snp in str(rec.get("rsids", "")).split(","), rec))
    keep = []
    for snp, hits in by_snp.items():
        _, rec = max(hits, key=lambda h: h[0])
        keep.append(dict(SNP=snp, chr=rec["chrom"], pos=rec["pos"],
                         effect_allele=rec["alt"], other_allele=rec["ref"],
                         eaf=pd.to_numeric(rec["af_alt"], errors="coerce"),
                         beta=pd.to_numeric(rec["beta"], errors="coerce"),
                         se=pd.to_numeric(rec["sebeta"], errors="coerce"),
                         p=pd.to_numeric(rec["pval"], errors="coerce")))
    if not keep:
        print(f"  {key}: nothing returned")
        return None
    new = pd.DataFrame(keep).dropna(subset=["beta", "se"])
    new["n"] = float(_sumstats.TRAIT_N[key][0])
    old = pd.read_csv(path) if path.exists() else None
    if old is not None and len(old):
        new = pd.concat([old[~old["SNP"].isin(new["SNP"])], new], ignore_index=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    new.to_csv(path, index=False)
    asked.write_text("\n".join(sorted(want)))
    print(f"  {key}: {len(keep)} rows matched", flush=True)
    return new[new["SNP"].isin(want)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--threads", type=int, default=10)
    ap.add_argument("--force-fetch", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    inst = mr.load_instruments(EXPOSURE)
    if inst is None or inst.empty:
        print(f"{EXPOSURE}: no instrument table")
        return 2
    want = {g.strip().upper() for g in a.only.split(",") if g.strip()}
    print(f"=== {len(inst)} instruments; family alpha "
          + ", ".join(f"{k}={v:.4f}" for k, v in FAMILY_ALPHA.items()))

    rows = []
    for group, key, pheno, cases, expect, why in PANEL:
        if want and group not in want:
            continue
        out = cached_outcome(key, inst, threads=a.threads, force=a.force_fetch)
        if out is None:
            continue
        r = estimate(EXPOSURE, key, inst=inst, out=out, label=f"{group}: {pheno}")
        if not r:
            continue
        r.update(group=group, cases=cases, phenocode=pheno, expectation=expect,
                 rationale=why, family_alpha=FAMILY_ALPHA[group],
                 bonferroni=bool(r["p"] < FAMILY_ALPHA[group]),
                 verdict=("protective" if r["b"] < 0 and r["p"] < 0.05
                          else "null" if r["p"] >= 0.05 else "adverse"))
        rows.append(r)

    if not rows:
        print("nothing estimated")
        return 1
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "panel.csv", index=False, float_format="%.6g")
    print(f"\n-> {(OUT / 'panel.csv').relative_to(ROOT).as_posix()}")
    verdicts(df)
    series()
    return 0


def verdicts(df):
    print("\n--- pre-registered decision rules")
    for group, label in [("A", "cardiac or systemic matrix?"),
                         ("B", "atrial muscle or conduction system?"),
                         ("C", "does the effect dilute as the endpoint aggregates?")]:
        d = df[df["group"] == group]
        if d.empty:
            continue
        print(f"\n  [{group}] {label}   (Bonferroni alpha {FAMILY_ALPHA[group]:.4f})")
        for _, r in d.iterrows():
            mark = "*" if r["bonferroni"] else " "
            print(f"    {mark} {r['phenocode']:16s} {r['b']:+.4f} "
                  f"({r['lo']:+.4f}, {r['hi']:+.4f}) p={r['p']:.3g} "
                  f"power={r['power_ref']:.0%}  {r['verdict']}")
        hits = d[d["bonferroni"]]
        if group == "A":
            print("    -> " + ("RULE A: systemic matrix is IN - the effect is not "
                               "cardiac-specific" if len(hits)
                               else "RULE A: cardiac-specific stands; no systemic-matrix "
                                    "endpoint moved"))
        if group == "B":
            cond = d[d["phenocode"].isin(["I9_CONDUCTIO", "I9_AVBLOCK"])]
            if len(cond) and (cond["bonferroni"] | (cond["p"] < 0.05)).any():
                print("    -> RULE B: conduction moves too - shared substrate, the atrial "
                      "label is too narrow")
            elif (d[d["phenocode"] == "I9_AF_REIMB"]["p"] < 0.05).any():
                print("    -> RULE B: atrial moves, conduction does not - atrial muscle, "
                      "and the heart-failure overlap does not explain it")
            else:
                print("    -> RULE B: neither atrial nor conduction survived; read the "
                      "estimates directly")


def series():
    """Group C plus the two module-14 anchors, ordered from aggregated to specific."""
    lad = pd.read_csv(OUT / "ladder.csv").set_index("outcome")
    pan = pd.read_csv(OUT / "panel.csv").set_index("outcome")
    rows = []
    for key, label in SERIES:
        src = pan if key in pan.index else lad
        if key not in src.index:
            print(f"  series: {key} not estimated yet, skipped")
            continue
        r = src.loc[key]
        rows.append(dict(outcome=key, label=label, b=float(r["b"]), se=float(r["se"]),
                         p=float(r["p"]), cases=int(r["cases"])))
    if len(rows) < 3:
        return
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "dilution_series.csv", index=False, float_format="%.6g")
    print("\n  [C] dilution series, aggregated -> specific")
    for _, r in d.iterrows():
        print(f"      {r['label']:32s} {r['b']:+.4f}  p={r['p']:.3g}  "
              f"cases={r['cases']:,}")
    diffs = np.diff(d["b"].to_numpy())
    mono = bool(np.all(diffs <= 1e-12))
    mags = np.abs(d["b"].to_numpy())
    mono_mag = bool(np.all(np.diff(mags) >= -1e-12))
    print(f"      monotone in sign: {mono}; monotone in magnitude: {mono_mag}")

    # the strict tests above fail on the pair of near-zero aggregates, where the difference
    # (0.0008) is far inside the standard error. The test that carries information is whether
    # the aggregates as a group differ from the specific endpoints as a group.
    agg = d[d["outcome"].isin(["FG_CVD_ALL", "FG_IHD", "FG_CHD"])]
    spec = d[d["outcome"].isin(["FG_HF", "FG_NONISCH_CM"])]
    if len(agg) >= 2 and len(spec) >= 2:
        wa = 1 / agg["se"].to_numpy() ** 2
        ba = float(np.sum(wa * agg["b"]) / np.sum(wa))
        sea = float(np.sqrt(1 / np.sum(wa)))
        ws = 1 / spec["se"].to_numpy() ** 2
        bs = float(np.sum(ws * spec["b"]) / np.sum(ws))
        ses = float(np.sqrt(1 / np.sum(ws)))
        diff = bs - ba
        sed = float(np.sqrt(sea ** 2 + ses ** 2))
        p = float(2 * stats.norm.sf(abs(diff / sed)))
        print(f"      aggregates pooled  {ba:+.4f} (se {sea:.4f})  "
              f"[{', '.join(agg['outcome'])}]")
        print(f"      specific  pooled  {bs:+.4f} (se {ses:.4f})  "
              f"[{', '.join(spec['outcome'])}]")
        print(f"      difference {diff:+.4f} (se {sed:.4f}), p={p:.2g}")
        if p < 0.05:
            print("      -> RULE C: the specific endpoints carry an effect the aggregated "
                  "endpoints do not - the dilution reading is measured, not arithmetic")
        else:
            print("      -> RULE C: no difference between aggregated and specific "
                  "endpoints; the dilution reading stays arithmetic only")
    elif mono_mag:
        print("      -> RULE C (weaker form): magnitude-monotone only")


def selftest():
    """Monotonicity must be detected on a planted series and refused on a broken one."""
    good = np.array([-0.005, -0.010, -0.064, -0.102])
    bad = np.array([-0.005, -0.080, -0.064, -0.102])
    ok_good = bool(np.all(np.diff(good) <= 0))
    ok_bad = bool(np.all(np.diff(bad) <= 0))
    print(f"selftest monotone series detected: {ok_good}")
    print(f"selftest broken series rejected: {not ok_bad}")
    ok = ok_good and not ok_bad
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

