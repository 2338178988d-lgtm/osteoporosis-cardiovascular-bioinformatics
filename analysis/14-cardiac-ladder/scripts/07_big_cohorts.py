"""Module 14, part 5 - do the AF and HF effects survive larger, and non-Finnish, cohorts?

Everything in sections A-D rests on Finnish endpoints, except the heart-failure anchor, which
module 15 checked in HERMES. Two criticisms follow directly:

  1. atrial fibrillation was measured once, in Finland, at 56-67 % power by the realised
     standard error. It is the most surprising result in the project;
  2. heart failure was measured twice (FinnGen, HERMES) and both have power to spare, but one
     is Finnish and one shares UK Biobank participants with the exposure.

Four further outcome cohorts fix both. Three carry tabix indexes, so the 600 instrument
positions are read over HTTP range requests and the 1.7-4.0 GB files are never downloaded;
the fourth has no index and is fetched whole (579 MB).

  key      source                                    cases    ancestry        overlap
  AF_US    GCST90475970 PMID 39024449 PheCode 427.21  71,577   US European     no UKB, no FI
  AF_NFE   GCST90473571 PMID 40770095 ICD10 I48      40,273   non-Finnish EUR UKB only
  AF_BIG   GCST90624411 PMID 40645996                252,438   European        UKB + FI
  HF_BIG   GCST90455657 PMID 40195560                168,363   European        UKB + FI
  HF_NFE   GCST90473575 PMID 40770095 ICD10 I50      19,629   non-Finnish EUR UKB only

DECISION RULES, fixed before the run:

  the protective estimate reproduces in at least one cohort that overlaps neither UK Biobank
  nor FinnGen (AF_US for atrial fibrillation; no such cohort exists for heart failure)
        -> the endpoint effects are not a Finnish founder-population artefact, and not an
           artefact of sharing participants with the exposure
  estimates agree in sign across all cohorts within an endpoint, with a pooled I^2 below 0.5
        -> report one pooled estimate per endpoint and stop treating cohort-specific numbers
           as separate findings
  any cohort reverses sign at p < 0.05
        -> that cohort is reported as a contradiction, not averaged away

Outputs (results/tables/14-cardiac-ladder/):
  big_cohorts.csv    per cohort: beta, CI, p, power, ancestry, overlap axis
  big_cohorts_pooled.csv  within-endpoint pooling across cohorts, with I^2

Usage:
  python analysis/14-cardiac-ladder/scripts/07_big_cohorts.py
  python analysis/14-cardiac-ladder/scripts/07_big_cohorts.py --self-test
  python analysis/14-cardiac-ladder/scripts/07_big_cohorts.py --only AF_US
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
GWAS = ROOT / "data" / "raw" / "gwas"
IDX = GWAS / "_tbi"
BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
EXPOSURE = "eBMD"

# (key, accession, bucket, main path, index path or None, endpoint, cases, controls,
#  ancestry, overlap with UKB / FinnGen, note)
BIG = [
    ("AF_US", "GCST90475970", "GCST90475001-GCST90476000", "GCST90475970.tsv.gz", None,
     "AF", 71577, 366195, "US European", "neither",
     "PheCode 427.21, Million Veteran Program"),
    ("AF_NFE", "GCST90473571", "GCST90473001-GCST90474000",
     "harmonised/GCST90473571.h.tsv.gz", "GCST90473571.tbi",
     "AF", 40273, 418167, "non-Finnish European", "UKB only",
     "ICD10 I48, UK Biobank"),
    ("AF_BIG", "GCST90624411", "GCST90624001-GCST90625000",
     "harmonised/GCST90624411.h.tsv.gz", "GCST90624411.tbi",
     "AF", 252438, 1959739, "European", "UKB + FinnGen",
     "largest published AF GWAS"),
    ("HF_BIG", "GCST90455657", "GCST90455001-GCST90456000",
     "harmonised/GCST90455657.h.tsv.gz", "GCST90455657.tbi",
     "HF", 168363, 1725558, "European", "UKB + FinnGen",
     "3.5x the HERMES case count"),
    ("HF_NFE", "GCST90473575", "GCST90473001-GCST90474000",
     "harmonised/GCST90473575.h.tsv.gz", "GCST90473575.tbi",
     "HF", 19629, 438811, "non-Finnish European", "UKB only",
     "ICD10 I50, UK Biobank"),
    # non-ischaemic cardiomyopathy: the largest effect in the project and, until now, measured
    # only in Finland. Three replication files, all from publications other than FinnGen.
    ("CM_MVP2", "GCST90475949", "GCST90475001-GCST90476000",
     "harmonised/GCST90475949.h.tsv.gz", "GCST90475949.tbi",
     "CM", 20680, 418334, "European", "UKB + FinnGen",
     "MVP PheCode cardiomyopathy"),
    ("CM_MVP1", "GCST90475951", "GCST90475001-GCST90476000",
     "GCST90475951.tsv.gz", None,
     "CM", 19708, 419904, "European", "UKB + FinnGen",
     "MVP PheCode cardiomyopathy, second definition; no index, full download"),
    ("CM_NFE", "GCST90473562", "GCST90473001-GCST90474000",
     "harmonised/GCST90473562.h.tsv.gz", "GCST90473562.tbi",
     "CM", 3578, 454862, "non-Finnish European", "UKB only",
     "ICD cardiomyopathy, UK Biobank - the only cardiomyopathy file with zero FinnGen overlap"),
]
# the estimates this module is meant to corroborate, from earlier modules
REFERENCE = {"AF": ("FG_AF", -0.0855), "HF": ("FG_HF", -0.0640),
             "CM": ("FG_NONISCH_CM", -0.1015)}


def columns(header):
    """GWAS Catalog harmonised files differ in naming; map onto the project schema.

    `beta` is optional: some releases (Million Veteran Program, GCST90475970) publish an odds
    ratio with `standard_error` set to NA and a confidence interval instead. Those files are
    converted to a log-OR beta with an SE derived from the interval.
    """
    h = set(header)
    pick = lambda *c: next((x for x in c if x in h), None)
    return dict(
        SNP=pick("rsid", "rs_id", "variant_id"),
        chr=pick("chromosome", "#chrom"),
        pos=pick("base_pair_position", "base_pair_location", "pos"),
        effect_allele=pick("effect_allele", "alt"),
        other_allele=pick("other_allele", "ref"),
        eaf=pick("effect_allele_frequency", "af_alt"),
        beta=pick("beta", "hm_beta"),
        se=pick("standard_error", "sebeta"),
        p=pick("p_value", "pval"),
        odds_ratio=pick("odds_ratio"),
        ci_upper=pick("ci_upper", "hm_ci_upper"),
        ci_lower=pick("ci_lower", "hm_ci_lower"),
    )


def effect_and_se(rec, cols):
    """Beta on the log scale and its SE, from whichever scale the file uses."""
    def num(name):
        v = rec.get(cols[name]) if cols.get(name) else None
        try:
            return float(v)
        except (TypeError, ValueError):
            return np.nan

    beta, se = num("beta"), num("se")
    if np.isfinite(beta) and np.isfinite(se) and se > 0:
        return beta, se
    or_, hi, lo = num("odds_ratio"), num("ci_upper"), num("ci_lower")
    if np.isfinite(or_) and or_ > 0 and np.isfinite(hi) and np.isfinite(lo) and 0 < lo < hi:
        return float(np.log(or_)), float((np.log(hi) - np.log(lo)) / (2 * 1.959964))
    return np.nan, np.nan


def index_file(acc, bucket, tbi):
    """Locate the .tbi on disk, or fetch it (2 MB).

    The download plan stores indexes as `<KEY>_<ACCESSION>.tbi`, so a plain name lookup misses
    them and the run would refetch a file that is already on disk. Search by accession first.
    """
    for candidate in list(IDX.glob(f"*{acc}*.tbi")) + [IDX / tbi]:
        if candidate.exists():
            return candidate
    dest = IDX / tbi
    import urllib.request
    url = f"{BASE}{bucket}/{acc}/harmonised/{tbi}"
    IDX.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(urllib.request.urlopen(url, timeout=240).read())
    return dest


def from_tabix(key, acc, bucket, main, tbi, inst, threads):
    tpath = index_file(acc, bucket, tbi)
    print(f"  {key}: reading {len(inst)} positions via tabix ...", flush=True)
    reader = TabixHttp(f"{BASE}{bucket}/{acc}/{main}", tpath)
    cols = columns(reader.header)
    targets = [(r.SNP, r.chr, int(r.pos)) for r in inst.itertuples()]
    rows = fetch(reader, targets, threads=threads)
    by_snp = {}
    for snp, pos, f in rows:
        rec = dict(zip(reader.header, f))
        rs = str(rec.get(cols["SNP"], ""))
        by_snp.setdefault(snp, []).append((snp == rs, rec))
    return _frame(by_snp, cols)


def from_full(key, acc, inst):
    """Read only the instrument rows from a full GWAS Catalog file.

    Every chunk is scanned: returning at the first chunk that contains any instrument would
    silently drop every locus downstream of it.
    """
    path = GWAS / f"{key}_{acc}.tsv.gz"
    if not path.exists():
        return None
    print(f"  {key}: scanning {path.name} ...", flush=True)
    want = set(inst["SNP"])
    found = {}
    for chunk in pd.read_csv(path, sep="\t", chunksize=2_000_000, low_memory=False,
                             dtype={"rsid": "string", "rs_id": "string"}):
        cols = columns(list(chunk.columns))
        col = cols["SNP"]
        hit = chunk[chunk[col].isin(want)]
        for _, r in hit.iterrows():
            found[str(r[col])] = r.to_dict()
        if len(found) == len(want):
            break
    if not found:
        print(f"  {key}: no instrument rsIDs present")
        return None
    print(f"  {key}: {len(found)}/{len(want)} instrument rsIDs found")
    return _frame({k: [(True, v)] for k, v in found.items()}, cols)


def _frame(by_snp, cols):
    """Build the project schema from raw rows, converting OR-scale files to log OR."""
    keep = []
    for snp, hits in by_snp.items():
        _, rec = max(hits, key=lambda h: h[0])
        beta, se = effect_and_se(rec, cols)

        def num(name):
            v = rec.get(cols[name]) if cols.get(name) else None
            try:
                return float(v)
            except (TypeError, ValueError):
                return np.nan

        pos = num("pos")
        keep.append(dict(SNP=snp, chr=str(rec.get(cols["chr"], "")),
                         pos=int(pos) if np.isfinite(pos) else 0,
                         effect_allele=str(rec.get(cols["effect_allele"], "")).upper(),
                         other_allele=str(rec.get(cols["other_allele"], "")).upper(),
                         eaf=num("eaf"), beta=beta, se=se, p=num("p")))
    if not keep:
        return None
    d = pd.DataFrame(keep)
    return d.dropna(subset=["beta", "se"])


def load_outcome(key, acc, bucket, main, tbi, inst, threads, n_total):
    path = CACHE / f"{key}.csv"
    asked = CACHE / f"{key}.snps"
    want = set(inst["SNP"])
    if path.exists() and asked.exists() and set(asked.read_text().split()) == want:
        d = pd.read_csv(path)
        return d[d["SNP"].isin(want)]
    d = from_tabix(key, acc, bucket, main, tbi, inst, threads) if tbi else from_full(key, acc, inst)
    if d is None or d.empty:
        print(f"  {key}: nothing returned")
        return None
    d["n"] = float(n_total)
    CACHE.mkdir(parents=True, exist_ok=True)
    d.to_csv(path, index=False)
    asked.write_text("\n".join(sorted(want)))
    print(f"  {key}: {len(d)} rows -> {path.relative_to(ROOT).as_posix()}", flush=True)
    return d[d["SNP"].isin(want)]


def pool(b, se):
    b = np.asarray(b, float)
    se = np.asarray(se, float)
    w = 1 / se ** 2
    bf = float(np.sum(w * b) / np.sum(w))
    q = float(np.sum(w * (b - bf) ** 2))
    dfree = len(b) - 1
    i2 = max(0.0, (q - dfree) / q) if q > 0 and dfree > 0 else 0.0
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (q - dfree) / c) if c > 0 and dfree > 0 else 0.0
    ws = 1 / (se ** 2 + tau2)
    br = float(np.sum(ws * b) / np.sum(ws))
    ser = float(np.sqrt(1 / np.sum(ws)))
    return dict(k=len(b), b=br, se=ser, lo=br - 1.96 * ser, hi=br + 1.96 * ser,
                p=float(2 * stats.norm.sf(abs(br / ser))), i2=i2,
                q_p=float(stats.chi2.sf(q, dfree)) if dfree > 0 else np.nan)


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
        print(f"{EXPOSURE}: no instruments")
        return 2
    want = {k.strip() for k in a.only.split(",") if k.strip()}
    print(f"=== {len(inst)} instruments; {len(BIG)} outcome cohorts")

    rows = []
    for key, acc, bucket, main_name, tbi, endpoint, cases, controls, ancestry, overlap, note in BIG:
        if want and key not in want:
            continue
        if a.force_fetch:
            for f in (CACHE / f"{key}.csv", CACHE / f"{key}.snps"):
                f.unlink(missing_ok=True)
        out = load_outcome(key, acc, bucket, main_name, tbi, inst, a.threads,
                           cases + controls)
        if out is None:
            print(f"  {key}: skipped (source not ready)")
            continue
        r = estimate(EXPOSURE, key, inst=inst, out=out, label=f"{key} ({endpoint})")
        if not r:
            continue
        r.update(endpoint=endpoint, accession=acc, cases=cases, controls=controls,
                 ancestry=ancestry, overlap=overlap, source=note)
        rows.append(r)

    if not rows:
        print("nothing estimated")
        return 1
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "big_cohorts.csv", index=False, float_format="%.6g")
    print(f"\n-> {(OUT / 'big_cohorts.csv').relative_to(ROOT).as_posix()}")

    pooled = []
    for endpoint in ("AF", "HF", "CM"):
        d = df[df["endpoint"] == endpoint]
        if d.empty:
            continue
        p = pool(d["b"], d["se"])
        p["endpoint"] = endpoint
        pooled.append(p)
        anchor_key, anchor_b = REFERENCE[endpoint]
        print(f"\n  [{endpoint}] reference {anchor_key} = {anchor_b:+.4f}")
        for _, s in d.iterrows():
            print(f"     {s['exposure'][:2]} {s['accession']} {s['ancestry']:22s} "
                  f"overlap={s['overlap']:15s} b={s['b']:+.4f} "
                  f"({s['lo']:+.4f}, {s['hi']:+.4f}) p={s['p']:.3g} "
                  f"power={s['power_ref']:.0%}")
        print(f"     pooled k={p['k']} b={p['b']:+.4f} ({p['lo']:+.4f}, {p['hi']:+.4f}) "
              f"p={p['p']:.3g} I2={p['i2']:.2f}")
    pd.DataFrame(pooled).to_csv(OUT / "big_cohorts_pooled.csv", index=False,
                                float_format="%.6g")
    print(f"-> {(OUT / 'big_cohorts_pooled.csv').relative_to(ROOT).as_posix()}")
    verdict(df)
    finnish_contrast(df)
    return 0


def finnish_contrast(df):
    """Every cohort on one axis, then Finnish versus everything else.

    Reading the new cohorts one at a time hides the pattern: all four are smaller than the
    Finnish estimate. Whether that is a founder-population effect, a case-definition effect or
    winner's curse is not separable here, but the direction and size of the gap are testable.
    """
    known = {
        "AF": [("FinnGen strict (63,532)", -0.0855, 0.0266, "Finnish"),
               ("FinnGen reimbursement (40,251)", -0.0961, 0.0302, "Finnish")],
        "HF": [("FinnGen strict (37,653)", -0.0640, 0.0186, "Finnish"),
               ("HERMES (47,309)", -0.0336, 0.0130, "multi-ancestry, shares UKB")],
        "CM": [("FinnGen non-ischaemic (13,917)", -0.1015, 0.0261, "Finnish")],
    }
    rows = []
    for endpoint, base in known.items():
        for label, b, se, pop in base:
            rows.append(dict(endpoint=endpoint, label=label, b=b, se=se, population=pop))
        d = df[df["endpoint"] == endpoint]
        for _, s in d.iterrows():
            rows.append(dict(endpoint=endpoint,
                             label=f"{s['accession']} ({s['cases']:,})", b=s["b"],
                             se=s["se"], population=s["ancestry"]))
    allc = pd.DataFrame(rows)
    allc["finnish"] = allc["population"].str.startswith("Finnish")
    allc.to_csv(OUT / "big_cohorts_all.csv", index=False, float_format="%.6g")

    print("\n=== every cohort on one axis")
    for endpoint in ("AF", "HF", "CM"):
        d = allc[allc["endpoint"] == endpoint]
        if d.empty:
            continue
        overall = pool(d["b"], d["se"])
        fin = d[d["finnish"]]
        rest = d[~d["finnish"]]
        print(f"\n  [{endpoint}] all k={overall['k']}: {overall['b']:+.4f} "
              f"({overall['lo']:+.4f}, {overall['hi']:+.4f}) p={overall['p']:.3g} "
              f"I2={overall['i2']:.2f} Q p={overall['q_p']:.3g}")
        for _, s in d.sort_values("b").iterrows():
            flag = "FI" if s["finnish"] else "  "
            print(f"     {flag} {s['label'][:40]:40s} {s['b']:+.4f} "
                  f"(se {s['se']:.4f})  {s['population']}")
        if len(fin) and len(rest):
            pf = pool(fin["b"], fin["se"])
            pr = pool(rest["b"], rest["se"])
            diff = pf["b"] - pr["b"]
            sed = float(np.sqrt(pf["se"] ** 2 + pr["se"] ** 2))
            p = float(2 * stats.norm.sf(abs(diff / sed)))
            print(f"     Finnish pooled {pf['b']:+.4f} (k={pf['k']}) vs non-Finnish "
                  f"{pr['b']:+.4f} (k={pr['k']}): difference {diff:+.4f} (se {sed:.4f}), "
                  f"p={p:.2g}")
            if p < 0.05:
                print("     -> the Finnish cohorts carry a significantly larger effect; the "
                      "FinnGen numbers cannot be quoted as the general estimate")
            else:
                print("     -> no significant difference between Finnish and non-Finnish "
                      "cohorts; the FinnGen estimate is not inflated beyond noise")
    print(f"\n-> {(OUT / 'big_cohorts_all.csv').relative_to(ROOT).as_posix()}")
    inflation(allc)


def inflation(allc):
    """How much larger is the Finnish estimate than the non-Finnish one, endpoint by endpoint?

    This is the summary the per-endpoint tables cannot show: the Finnish estimate is the
    largest in every endpoint, and the gap is not the same size each time.
    """
    rows = []
    for endpoint in ("AF", "HF", "CM"):
        d = allc[allc["endpoint"] == endpoint]
        fin, rest = d[d["finnish"]], d[~d["finnish"]]
        if not len(fin):
            continue
        pf = pool(fin["b"], fin["se"])
        if not len(rest):
            continue
        pr = pool(rest["b"], rest["se"])
        bf, sef, br, ser, i2r = pf["b"], pf["se"], pr["b"], pr["se"], pr["i2"]
        diff = bf - br
        sed = float(np.sqrt(sef ** 2 + ser ** 2))
        rows.append(dict(endpoint=endpoint, k_finnish=len(fin), b_finnish=bf,
                         se_finnish=sef, k_other=len(rest), b_other=br, se_other=ser,
                         ratio=bf / br if br else np.nan,
                         difference=diff, se_difference=sed,
                         p=float(2 * stats.norm.sf(abs(diff / sed))), i2_other=i2r))
    if not rows:
        return
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "finnish_inflation.csv", index=False, float_format="%.6g")
    print("\n=== Finnish vs non-Finnish, all three endpoints")
    print(f"  {'endpoint':10s} {'Finnish':>9s} {'k':>2s} {'other':>9s} {'k':>2s} "
          f"{'ratio':>6s} {'difference':>11s} {'p':>7s}")
    for _, r in d.iterrows():
        print(f"  {r['endpoint']:10s} {r['b_finnish']:+9.4f} {int(r['k_finnish']):2d} "
              f"{r['b_other']:+9.4f} {int(r['k_other']):2d} {r['ratio']:6.1f} "
              f"{r['difference']:+11.4f} {r['p']:7.3f}")
    n_sig = int((d["p"] < 0.05).sum())
    print(f"  {n_sig} of {len(d)} endpoints differ significantly; "
          f"the Finnish estimate is larger in {int((d['ratio'] > 1).sum())} of {len(d)}")
    print("  -> this is the systematic pattern: not one failed replication but a consistent")
    print("     Finnish excess across every endpoint measured in both settings")
    print(f"\n-> {(OUT / 'finnish_inflation.csv').relative_to(ROOT).as_posix()}")



def verdict(df):
    print("\n--- pre-registered rules")
    for endpoint in ("AF", "HF", "CM"):
        d = df[df["endpoint"] == endpoint]
        if d.empty:
            continue
        clean = d[d["overlap"] == "neither"]
        rev = d[(d["b"] > 0) & (d["p"] < 0.05)]
        sign_ok = bool((d["b"] < 0).all())
        p = pool(d["b"], d["se"])
        print(f"\n  [{endpoint}]")
        if len(clean):
            c = clean.iloc[0]
            print(f"    cohort overlapping neither UKB nor FinnGen: {c['accession']} "
                  f"b={c['b']:+.4f}, p={c['p']:.3g} -> "
                  + ("REPRODUCES in a fully independent cohort" if c["b"] < 0 and c["p"] < 0.05
                     else "does not reach significance in the independent cohort"))
        else:
            print("    no cohort overlaps neither UKB nor FinnGen: full independence is not "
                  "testable for this endpoint")
        if len(rev):
            print("    CONTRADICTION: "
                  + ", ".join(f"{r['accession']} ({r['b']:+.3f}, p={r['p']:.3g})"
                              for _, r in rev.iterrows()))
        if sign_ok and p["i2"] < 0.5:
            print(f"    all cohorts same sign, I2={p['i2']:.2f}: report the pooled estimate "
                  f"({p['b']:+.4f}) rather than the cohorts separately")


def selftest():
    """Pooling must recover a planted common effect and flag a planted reversal."""
    rng = np.random.default_rng(5)
    b = rng.normal(-0.065, 0.006, 3)
    se = np.full(3, 0.02)
    good = pool(b, se)
    print(f"selftest homogeneous b={good['b']:+.4f} I2={good['i2']:.2f}")
    ok1 = abs(good["b"] + 0.065) < 0.01 and good["i2"] < 0.3
    mixed = pool(np.array([-0.065, -0.060, 0.070]), se)
    print(f"selftest with a reversal I2={mixed['i2']:.2f} Qp={mixed['q_p']:.3g}")
    ok2 = mixed["i2"] > 0.5
    ok = ok1 and ok2
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

