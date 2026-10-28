"""Module 14, part 3 鈥?three checks a reviewer will ask for before believing the panel.

None of these needs new data.

  1. SHARED-LOCUS REMOVAL
     The twelve panel endpoints and the three protective module-14 rungs all come from the
     same 600 instruments, and heart failure is a highly polygenic trait whose own top hits
     could coincide with bone-density loci. If an instrument also associates with heart
     failure, its panel estimate partly reflects heart failure, not the panel endpoint.
     Module 05 did this for the heart-failure anchor; the panel has to do the same. Any eBMD
     instrument within 500 kb of a genome-wide-significant heart-failure variant is removed
     and every endpoint is re-estimated.

  2. STEIGER FILTERING
     An instrument that explains more variance in the outcome than in the exposure points the
     wrong way. Module 05 filtered for heart failure; the panel is unfiltered. For each
     endpoint the per-variant R^2 is computed on both sides and variants with
     R^2_outcome > R^2_exposure are dropped.

  3. EFFECTIVE NUMBER OF TESTS
     The panel endpoints are correlated (per-variant Wald ratios correlate up to 0.69), so
     Bonferroni's 0.0125 is conservative. The effective number of independent tests is
     estimated from the eigenvalue spectrum of the estimate correlation matrix (Li & Ji 2005)
     and a corrected threshold is reported next to the Bonferroni one.

All three are recomputed from the caches in data/processed/mr/outcome_snps/, so the run is
offline.

Outputs (results/tables/14-cardiac-ladder/):
  panel_checks.csv      per endpoint: base, shared-locus-removed, Steiger-filtered estimates
  effective_tests.csv   Meff and the corrected threshold for the panel and the ladder

Usage:
  python analysis/14-cardiac-ladder/scripts/05_panel_checks.py
  python analysis/14-cardiac-ladder/scripts/05_panel_checks.py --self-test
"""
import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "analysis" / "scripts"))
sys.path.insert(0, str(ROOT / "analysis" / "13-robustness" / "scripts"))

from _common import mr  # noqa: E402
import _sumstats  # noqa: E402
from _tabix_http import TabixHttp, fetch  # noqa: E402

# module 05's Steiger implementation, so the definition of "exposure-dominant" does not
# drift between modules
import importlib.util as _ilu  # noqa: E402
_spec = _ilu.spec_from_file_location(
    "_hf_followup", ROOT / "analysis" / "05-hf-followup" / "scripts" / "01_hf_followup.py")
_hf = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_hf)
steiger_filter = _hf.steiger_filter

OUT = ROOT / "results" / "tables" / "14-cardiac-ladder"
CACHE = ROOT / "data" / "processed" / "mr" / "outcome_snps"
FINNGEN_DIR = ROOT / "data" / "raw" / "finngen"
FG_URL = ("https://storage.googleapis.com/finngen-public-data-r12/summary_stats/"
          "release/finngen_R12_{}.gz")
EXPOSURE = "eBMD"
WINDOW_KB = 500

PANEL = ["FG_VARICVE", "FG_VTE", "FG_PULMEMB", "FG_AORTANEUR", "FG_OTHARR",
         "FG_AF_REIMB", "FG_PAROXTAC", "FG_CONDUCTIO", "FG_AVBLOCK",
         "FG_CVD_ALL", "FG_IHD", "FG_CHD"]
LADDER_KEYS = ["FG_HYPERTENS", "FG_MI_STRICT", "FG_AF", "FG_VHD", "FG_CAVS",
               "FG_HYPERTENS_HD", "FG_NONISCH_CM", "FG_CARDMYO"]


def load_from_cache(key):
    path = CACHE / f"{key}.csv"
    if not path.exists():
        return None
    d = pd.read_csv(path)
    return d if len(d) else None


def hf_top_hits():
    """Genome-wide-significant heart-failure variants, from the local FinnGen file.

    Read from `data/raw/finngen/finngen_R12_I9_HEARTFAIL.gz`, which module 05 already
    downloaded and verified. An earlier version of this function queried the remote tabix
    index instead; that index was never fetched for this endpoint (only "strict" HF was), so
    the query silently returned nothing and every "-N loci" in the output was actually zero.
    """
    path = CACHE / "FG_HF_TOPHITS.csv"
    if path.exists():
        return pd.read_csv(path)
    src = FINNGEN_DIR / "finngen_R12_I9_HEARTFAIL.gz"
    if not src.exists():
        print("  heart-failure file absent; shared-locus removal skipped")
        return None
    rows = []
    with gzip.open(src, "rt") as f:
        head = f.readline().rstrip("\n").lstrip("#").split("\t")
        idx = {c: i for i, c in enumerate(head)}
        for line in f:
            p = line.split("\t")
            try:
                if float(p[idx["pval"]]) < 5e-8:
                    rows.append(dict(SNP=p[idx["rsids"]].split(",")[0],
                                     chr=p[idx["chrom"]], pos=int(p[idx["pos"]]),
                                     p=float(p[idx["pval"]])))
            except (ValueError, IndexError):
                continue
    if not rows:
        return None
    d = pd.DataFrame(rows).sort_values("p").drop_duplicates("SNP")
    d.to_csv(path, index=False)
    print(f"  heart-failure top hits: {len(d)}")
    return d


def drop_shared(inst, hits):
    """eBMD instruments within WINDOW_KB of a heart-failure top hit."""
    if hits is None or hits.empty:
        return inst, 0
    h = hits.copy()
    h["chr"] = h["chr"].astype(str).str.replace(r"\.0$", "", regex=True)
    chr_i = inst["chr"].astype(float).astype(int).astype(str)
    pos_i = inst["pos"].astype(float)
    keep = np.ones(len(inst), bool)
    by_chr = {c: g["pos"].to_numpy(float) for c, g in h.groupby("chr")}
    for i, (c, p) in enumerate(zip(chr_i, pos_i)):
        arr = by_chr.get(c)
        if arr is None or not len(arr):
            continue
        if np.min(np.abs(arr - p)) <= WINDOW_KB * 1000:
            keep[i] = False
    return inst[keep].reset_index(drop=True), int((~keep).sum())


def steiger_on_harmonised(d):
    """Steiger filter on an already-harmonised table, using module 05's estimator."""
    if len(d) < 12:
        return d, 0
    keep, _, _ = steiger_filter(d)
    return d[keep], int((~keep).sum())


def steiger(inst, out):
    """Keep variants that explain more variance in the exposure than in the outcome."""
    return steiger_on_harmonised(mr.harmonise(inst, out))


def ivw_from(d):
    if len(d) < 8:
        return None
    return mr.ivw(d["bx"].to_numpy(float), d["by"].to_numpy(float),
                  d["sx"].to_numpy(float), d["sy"].to_numpy(float), random=True)


def run(inst, key, label):
    out = load_from_cache(key)
    if out is None:
        print(f"  {key}: no cached rows, skipped")
        return None
    base = ivw_from(mr.harmonise(inst, out))
    trimmed, n_shared = drop_shared(inst, HITS)
    trim = ivw_from(mr.harmonise(trimmed, out))
    filt, n_steiger = steiger(inst, out)
    st = ivw_from(filt)
    row = dict(outcome=key, label=label,
               b_base=base["b"] if base else np.nan,
               p_base=base["p"] if base else np.nan,
               nsnp_base=len(mr.harmonise(inst, out)),
               n_shared_removed=n_shared, n_instruments_left=len(trimmed),
               b_trimmed=trim["b"] if trim else np.nan,
               p_trimmed=trim["p"] if trim else np.nan,
               n_steiger_removed=n_steiger,
               b_steiger=st["b"] if st else np.nan,
               p_steiger=st["p"] if st else np.nan,
               nsnp_steiger=len(filt) if st else 0)
    row["flip_under_trim"] = bool(
        np.isfinite(row["b_base"]) and np.isfinite(row["b_trimmed"])
        and np.sign(row["b_base"]) != np.sign(row["b_trimmed"]))
    print(f"  {label:22s} base {row['b_base']:+.4f} (p={row['p_base']:.3g})  "
          f"trimmed {row['b_trimmed']:+.4f} (p={row['p_trimmed']:.3g}, "
          f"-{n_shared} loci)  steiger {row['b_steiger']:+.4f} "
          f"(p={row['p_steiger']:.3g}, -{n_steiger} SNPs)")
    return row


def effective_tests(keys):
    """Li & Ji (2005) effective number of independent tests from the correlation matrix."""
    series = []
    for k in keys:
        out = load_from_cache(k)
        if out is None:
            continue
        d = mr.harmonise(INST, out).set_index("SNP")
        series.append((k, d["by"] / d["bx"]))
    if len(series) < 3:
        return None
    w = pd.concat([s for _, s in series], axis=1)
    w.columns = [k for k, _ in series]
    c = w.corr().to_numpy()
    c = np.nan_to_num(c, nan=0.0)
    np.fill_diagonal(c, 1.0)
    eig = np.sort(np.linalg.eigvalsh(c))[::-1]
    eig = np.clip(eig, 0, None)
    meff = float(1 + (np.sum(eig > 1) - 1) * 0 + np.sum((eig - 1) * (eig > 1)))
    meff = float(max(1.0, np.sum(eig > 1) + np.sum(eig[eig > 1] - 1)))
    return dict(k=len(series), meff=meff, alpha_bonf=0.05 / len(series),
                alpha_meff=0.05 / meff), c, [k for k, _ in series]


def main():
    global HITS, INST
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    INST = mr.load_instruments(EXPOSURE)
    if INST is None or INST.empty:
        print(f"{EXPOSURE}: no instruments")
        return 2
    print(f"=== {len(INST)} instruments")
    print("=== heart-failure top hits (for shared-locus removal)")
    HITS = hf_top_hits()

    rows = []
    print("\n=== per-endpoint checks")
    for key in PANEL:
        r = run(INST, key, key.replace("FG_", "").lower())
        if r:
            rows.append(r)
    for key in LADDER_KEYS:
        r = run(INST, key, key.replace("FG_", "").lower())
        if r:
            rows.append(r)
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(OUT / "panel_checks.csv", index=False, float_format="%.6g")
        print(f"\n-> {(OUT / 'panel_checks.csv').relative_to(ROOT).as_posix()}")
        n_flip = int(df["flip_under_trim"].sum())
        print(f"  estimates that change sign when heart-failure loci are removed: {n_flip}")

    print("\n=== effective number of tests")
    res = effective_tests(PANEL + LADDER_KEYS)
    if res:
        info, c, keys = res
        pd.DataFrame([info]).to_csv(OUT / "effective_tests.csv", index=False)
        print(f"  k={info['k']} endpoints, Meff={info['meff']:.2f}")
        print(f"  Bonferroni alpha {info['alpha_bonf']:.4f} -> "
              f"Meff-corrected {info['alpha_meff']:.4f}")
        print(f"  -> {(OUT / 'effective_tests.csv').relative_to(ROOT).as_posix()}")
        pd.DataFrame(c, index=keys, columns=keys).to_csv(
            OUT / "endpoint_correlation.csv", float_format="%.4g")
    return 0


def selftest():
    """The two filters must remove what they claim and leave a clean set alone."""
    rng = np.random.default_rng(3)
    n = 300
    inst = pd.DataFrame(dict(SNP=[f"rs{i}" for i in range(n)],
                             chr=np.repeat(np.arange(1, 4), n // 3)[:n].astype(float),
                             pos=rng.integers(1_000_000, 90_000_000, n).astype(float)))
    hits = pd.DataFrame(dict(SNP=["rs1"], chr=[1.0], pos=[inst.loc[1, "pos"]]))
    kept, removed = drop_shared(inst, hits)
    print(f"selftest shared-locus removal: removed {removed}, kept {len(kept)}")
    ok1 = removed >= 1 and len(kept) == n - removed

    # realistic scale: se ~ 1/sqrt(2f(1-f)n) with n = 400,000 and f ~ 0.3
    d = pd.DataFrame(dict(bx=rng.normal(0.05, 0.01, 200), sx=np.full(200, 0.0025),
                          n_x=np.full(200, 400000.0), by=rng.normal(0, 0.003, 200),
                          sy=np.full(200, 0.0025), n_y=np.full(200, 400000.0)))
    # 100 variants that explain far more variance in the outcome than in the exposure
    d.loc[:99, "by"] = rng.normal(0.08, 0.003, 100)
    _, n_removed = steiger_on_harmonised(d)
    ok2 = n_removed > 50
    print(f"selftest steiger: removed {n_removed} of 200 (100 planted violators, expect >50)")
    ok = ok1 and ok2
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

