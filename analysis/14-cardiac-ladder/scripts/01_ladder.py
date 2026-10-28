"""Module 14 - where on the cardiac-load ladder does the bone-density effect sit?

Every statement about mechanism in modules 02-13 is negative or inferred: not coronary
(the HF+CHD subtype is null, MVMR on CAD does not attenuate), not circulating (ten mediators
null), not body size, not any annotated pathway, and LV mass moves by only -0.041 SD. The
one positive localisation - "the myocardium" - was reached by elimination, never measured,
because FinnGen R12 has no "heart failure without coronary disease" endpoint.

That endpoint does not exist, but six others do, and each sits at a defined rung of the
path from vascular load to myocardial failure:

    hypertension (154,630) -> hypertensive heart disease (11,728) -> heart failure (37,653)
    atrial fibrillation (63,532)          non-ischaemic cardiomyopathy (13,917)
    valvular heart disease (28,602)       calcific aortic stenosis, operated (12,418)
    myocardial infarction, strict (31,666)                       cardiomyopathy (7,649)

Same exposure (UK Biobank heel-ultrasound eBMD, GCST90029004), same 600 instruments, same
estimator, and every outcome is Finnish - no UK Biobank overlap anywhere.

DECISION RULES, written before the run (the outcome is not chosen after seeing it):

  A. non-ischaemic cardiomyopathy AND hypertensive heart disease are both protective
     (p < 0.05, same sign as heart failure), while hypertension and myocardial infarction
     are null
        -> the myocardial route is measured, not inferred; the paper can say "bone-myocardium
           axis" and the mechanism section has positive content.
  B. only heart failure survives, and the well-powered non-coronary rungs (atrial
     fibrillation 63,532 cases; hypertension 154,630) are null
        -> heart-failure specific but not localisable; the mechanism section must be rewritten
           as "could not be localised" and the claim drops one rank. The atrial-fibrillation
           null at 99 % power is itself a strong statement (not a rhythm effect).
  C. calcific aortic stenosis is protective
        -> the ectopic-calcification / matrix-vesicle reading becomes a competitor to the
           myocardial one and must be argued in the discussion, not dismissed.

Multiplicity: eight endpoints tested, Bonferroni alpha = 0.05/8 = 0.00625, reported alongside
the unadjusted p. The heart-failure anchor (FG_HF, already local) is run in the same table so
the ladder is read against an estimate computed identically, but it is not counted.

Inputs: eight FinnGen R12 `.tbi` indexes in data/raw/finngen/. The 815 MB summary file for
each endpoint is never downloaded - variant rows are read over HTTP range requests
(analysis/scripts/_tabix_http.py) and cached in data/processed/mr/outcome_snps/, so a rerun
is offline.

Outputs (results/tables/14-cardiac-ladder/):
  ladder.csv      one row per endpoint: beta, CI, p, power, Bonferroni verdict
  ladder_loo.csv  leave-one-out range per endpoint

Usage:
  python analysis/14-cardiac-ladder/scripts/01_ladder.py
  python analysis/14-cardiac-ladder/scripts/01_ladder.py --self-test
  python analysis/14-cardiac-ladder/scripts/01_ladder.py --only FG_AF
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

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

# (key, phenocode, cases, rung, what the rung tests)
LADDER = [
    ("FG_HF", "I9_HEARTFAIL", 37653, "heart failure", "anchor - already local"),
    ("FG_HYPERTENS", "I9_HYPTENS", 154630, "load",
     "the vascular load itself; a null means the effect is not working through blood pressure"),
    ("FG_MI_STRICT", "I9_MI_STRICT", 31666, "coronary",
     "atherosclerotic endpoint, well powered"),
    ("FG_AF", "I9_AF", 63532, "rhythm", "large, non-coronary, non-ventricular"),
    ("FG_VHD", "I9_VHD", 28602, "valve", "structural heart disease, non-coronary"),
    ("FG_CAVS", "I9_CAVS_OPERATED", 12418, "valve/calcific",
     "ectopic calcification - the calcification-paradox competitor"),
    ("FG_HYPERTENS_HD", "I9_HYPTENSHD", 11728, "load -> myocardium", "the transition rung"),
    ("FG_NONISCH_CM", "I9_NONISCHCARDMYOP", 13917, "myocardium",
     "the usable stand-in for 'heart failure without coronary disease'"),
    ("FG_CARDMYO", "I9_CARDMYO", 7649, "myocardium", "sensitivity, all cardiomyopathy"),
]
LADDER_KEYS = [r[0] for r in LADDER if r[0] != "FG_HF"]
ALPHA = 0.05 / len(LADDER_KEYS)


def cached_outcome(key, inst, threads=8, force=False):
    """Instrument-SNP rows for a remote endpoint, cached on disk after the first fetch.

    A companion `.snps` file records which variants were asked for. Without it the cache
    could only be validated by demanding that it hold every requested variant, which never
    happens - some instrument positions are absent from the Finnish release - so every rerun
    would silently re-fetch all 600 positions.
    """
    path = CACHE / f"{key}.csv"
    asked = CACHE / f"{key}.snps"
    want = set(inst["SNP"])
    if path.exists() and not force:
        got = pd.read_csv(path)
        if asked.exists() and set(asked.read_text().split()) == want:
            return got[got["SNP"].isin(want)]
    pheno = _sumstats.FILES[key].replace("finngen_R12_", "").replace(".gz", "")
    tbi = FINNGEN_DIR / f"finngen_R12_{pheno}.gz.tbi"
    if not tbi.exists():
        print(f"  {key}: index {tbi.name} absent - run download_plan.py --category ladder")
        return None
    print(f"  {key}: fetching {len(inst)} positions via tabix ...", flush=True)
    reader = TabixHttp(FG_URL.format(pheno), tbi)
    cols = reader.header
    targets = [(r.SNP, r.chr, int(r.pos)) for r in inst.itertuples()]
    rows = fetch(reader, targets, threads=threads)
    by_snp = {}
    for snp, pos, f in rows:
        rec = dict(zip(cols, f))
        rsids = str(rec.get("rsids", "")).split(",")
        by_snp.setdefault(snp, []).append((snp in rsids, rec))
    keep = []
    for snp, hits in by_snp.items():
        exact, rec = max(hits, key=lambda h: h[0])
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
    print(f"  {key}: {len(keep)} rows matched -> {path.relative_to(ROOT).as_posix()}",
          flush=True)
    return new[new["SNP"].isin(want)]


def leave_one_out(inst, out):
    d = mr.harmonise(inst, out)
    if len(d) < 10:
        return None
    bx, by = d["bx"].to_numpy(float), d["by"].to_numpy(float)
    sx, sy = d["sx"].to_numpy(float), d["sy"].to_numpy(float)
    b = []
    for i in range(len(d)):
        m = np.ones(len(d), bool)
        m[i] = False
        b.append(mr.ivw(bx[m], by[m], sx[m], sy[m], random=True)["b"])
    return dict(nsnp=len(d), loo_min=float(np.min(b)), loo_max=float(np.max(b)),
                loo_crosses_zero=bool(np.min(b) < 0 < np.max(b)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--force-fetch", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    inst = mr.load_instruments(EXPOSURE)
    if inst is None or inst.empty:
        print(f"{EXPOSURE}: no instrument table - run 02_extract_instruments.py")
        return 2
    print(f"=== exposure {EXPOSURE}: {len(inst)} instruments, min F={inst['F'].min():.1f}")
    print(f"=== Bonferroni alpha over {len(LADDER_KEYS)} rungs: {ALPHA:.5f}")

    want = {k.strip() for k in a.only.split(",") if k.strip()}
    rows, loos = [], []
    for key, pheno, cases, rung, note in LADDER:
        if want and key not in want:
            continue
        out = None
        if not _sumstats.path_of(key).exists():
            out = cached_outcome(key, inst, threads=a.threads, force=a.force_fetch)
            if out is None:
                continue
        r = estimate(EXPOSURE, key, inst=inst, out=out, label=f"{rung}: {pheno}")
        if not r:
            continue
        r.update(rung=rung, cases=cases, phenocode=pheno, tests=note,
                 is_anchor=key == "FG_HF", bonferroni_alpha=ALPHA,
                 significant_bonferroni=bool(key != "FG_HF" and r["p"] < ALPHA))
        rows.append(r)
        lo = leave_one_out(inst, out if out is not None
                           else mr.outcome_at_snps(key, inst["SNP"].tolist()))
        if lo:
            loos.append(dict(outcome=key, **lo))

    if not rows:
        print("no rung estimated")
        return 1
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "ladder.csv", index=False, float_format="%.6g")
    print(f"\n-> {(OUT / 'ladder.csv').relative_to(ROOT).as_posix()}")
    if loos:
        pd.DataFrame(loos).to_csv(OUT / "ladder_loo.csv", index=False, float_format="%.6g")
        print(f"-> {(OUT / 'ladder_loo.csv').relative_to(ROOT).as_posix()}")
    verdict(df)
    return 0


def verdict(df):
    """Report each pre-registered rule against what was estimated, then the observed pattern.

    The rules were fixed before the run; this function does not reinterpret them after the
    fact. Where none of them fires, it says so and states the pattern separately, labelled as
    observed rather than pre-registered.
    """
    d = df[~df["is_anchor"]].set_index("outcome")
    print("\n--- pre-registered rules")
    need = {"FG_NONISCH_CM", "FG_HYPERTENS_HD", "FG_HYPERTENS", "FG_MI_STRICT"}
    if not need.issubset(d.index):
        print("  A/B not decidable: a required rung is missing")
    else:
        cm, hhd = d.loc["FG_NONISCH_CM"], d.loc["FG_HYPERTENS_HD"]
        htn, mi = d.loc["FG_HYPERTENS"], d.loc["FG_MI_STRICT"]
        a_myo = cm["b"] < 0 and cm["p"] < 0.05 and hhd["b"] < 0 and hhd["p"] < 0.05
        a_load = htn["p"] >= 0.05 and mi["p"] >= 0.05
        print(f"  RULE A  myocardial route measured: "
              f"{'MET' if (a_myo and a_load) else 'NOT MET'} "
              f"(NONISCH_CM b={cm['b']:+.4f} p={cm['p']:.3g}; "
              f"HYPERTENS_HD b={hhd['b']:+.4f} p={hhd['p']:.3g} "
              f"power={hhd['power_ref']:.0%}; load+coronary null: {a_load})")
        sig = d[d["significant_bonferroni"]]
        print(f"  RULE B  only heart failure survives: "
              f"{'MET' if len(sig) == 0 else 'NOT MET'} - {len(sig)} rung(s) survive "
              f"Bonferroni ({ALPHA:.4f})")
    if "FG_CAVS" in d.index:
        c = d.loc["FG_CAVS"]
        print(f"  RULE C  calcific aortic stenosis protective: "
              f"{'MET' if (c['b'] < 0 and c['p'] < 0.05) else 'NOT MET'} "
              f"(b={c['b']:+.4f} p={c['p']:.3g} power={c['power_ref']:.0%})")

    print("\n--- observed pattern (not pre-registered)")
    sig = d[d["significant_bonferroni"]].sort_values("b")
    if len(sig):
        print("  protective and Bonferroni-significant: "
              + ", ".join(f"{k} ({v:+.3f}, p={d.loc[k, 'p']:.2g}, "
                          f"power {d.loc[k, 'power_ref']:.0%})"
                          for k, v in sig["b"].items()))
    nulls = d[(~d["significant_bonferroni"])].sort_values("power_ref")
    print("  null: " + ", ".join(f"{k} ({v:+.3f}, p={d.loc[k, 'p']:.2g}, "
                                 f"power {d.loc[k, 'power_ref']:.0%})"
                                 for k, v in nulls["b"].items()))
    wellpowered_nulls = nulls[nulls["power_ref"] >= 0.65]
    if len(wellpowered_nulls):
        print("  the null rungs that carry enough power to mean something: "
              + ", ".join(wellpowered_nulls.index))


def selftest():
    """The estimator must recover a planted -0.064 on a 592-SNP set and leave a null alone."""
    rng = np.random.default_rng(11)
    n = 592
    bx = rng.normal(0.035, 0.02, n)
    sx = np.full(n, 0.0035)
    sy = np.full(n, 0.0125)
    res = []
    for tag, by in [("planted -0.064", -0.064 * bx + rng.normal(0, sy)),
                    ("null", rng.normal(0, sy))]:
        r = mr.ivw(bx, by, sx, sy, random=True)
        res.append((tag, r["b"], r["p"]))
        print(f"selftest {tag:16s} b={r['b']:+.4f} p={r['p']:.3g}")
    ok = (abs(res[0][1] + 0.064) < 0.02 and res[0][2] < 0.01 and res[1][2] > 0.01)
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

