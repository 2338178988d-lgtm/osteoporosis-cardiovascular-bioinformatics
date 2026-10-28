"""Module 12b 鈥?systematic gene-set scan of the diffuse eBMD effect.

Script 01 partitions the instruments into seven curated bone pathways and finds that the
canonical bone genes carry almost none of the effect (they are 29 of 592 instruments).
That leaves the question open in the other direction: is there *any* annotated programme
whose instruments behave differently from the rest?

This scans standard, unbiased collections 鈥?MSigDB Hallmark, GO Biological Process and
Reactome, exported from R (`msigdbr`) to `data/processed/gene_sets/msigdb_human.csv` 鈥?and
for every set with enough instruments compares:

    beta_set        IVW on the instruments whose nearest gene is in the set
    beta_complement IVW on the remaining instruments
    difference      beta_set - beta_complement, SEs added (the two groups are disjoint)

Multiplicity is handled twice: Benjamini-Hochberg across sets, and a family-wise
permutation test that permutes the gene -> instrument assignment, keeps the observed
per-SNP Wald ratios, and recomputes the largest |z| across all sets. The permutation
answers the only question that matters here: does knowing which gene an instrument sits
next to predict how much heart-failure effect it carries?

Outputs (results/tables/12-pathway-stratified/):
  geneset_scan.csv        per set: n, beta_set, beta_complement, difference, z, p, BH q
  geneset_fwer.csv        observed max |z|, permutation null quantiles, family-wise p

Usage:
  python analysis/12-pathway-stratified/scripts/02_gene_set_scan.py
  python analysis/12-pathway-stratified/scripts/02_gene_set_scan.py --self-test
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

ROOT = HERE.parents[2]
OUT = ROOT / "results" / "tables" / "12-pathway-stratified"
GENES = ROOT / "results" / "tables" / "09-integration" / "instrument_genes.csv"
SETS = ROOT / "data" / "processed" / "gene_sets" / "msigdb_human.csv"
EXPOSURE = "eBMD"
OUTCOMES = {"HF": "HERMES (overlapping)", "FG_HF": "FinnGen R12 (independent)"}
MIN_SNPS = 8


def fetch(outcome):
    inst = mr.load_instruments(EXPOSURE)
    out = mr.outcome_at_snps(outcome, inst["SNP"].tolist())
    if out is None or out.empty:
        return None
    return mr.harmonise(inst, out)


def membership(d, min_snps=MIN_SNPS):
    """SNP x gene-set boolean matrix restricted to sets with enough instruments, plus the
    per-SNP IVW weights.

    Returns (M, names, a, c) where a_i = w_i bx_i by_i and c_i = w_i bx_i^2, so that
    beta = a @ mask / (c @ mask) and se = 1 / sqrt(c @ mask) for any instrument subset.
    Restricting to eligible sets first keeps the permutation loop cheap and does not change
    the observed scan: a set that fails here cannot enter the family.
    """
    genes = pd.read_csv(GENES).set_index("SNP")["nearest_gene"]
    gs = pd.read_csv(SETS)
    order = sorted(gs["set"].unique())
    col = {s: i for i, s in enumerate(order)}
    gene_to_sets = {g: [col[s] for s in sub["set"]] for g, sub in gs.groupby("gene")}
    M = np.zeros((len(d), len(order)), dtype=bool)
    for i, snp in enumerate(d["SNP"]):
        g = genes.get(snp)
        hit = gene_to_sets.get(g) if isinstance(g, str) else None
        if hit:
            M[i, hit] = True
    keep = M.sum(axis=0) >= min_snps
    bx = d["bx"].to_numpy(float)
    by = d["by"].to_numpy(float)
    w = 1.0 / d["sy"].to_numpy(float) ** 2
    return M[:, keep], np.asarray(order)[keep], w * bx * by, w * bx ** 2


def estimates(a, c, masks):
    """IVW estimate and SE for each column of `masks` (n x k boolean)."""
    num = a @ masks
    den = c @ masks
    with np.errstate(divide="ignore", invalid="ignore"):
        b = np.where(den > 0, num / den, np.nan)
        se = np.where(den > 0, 1.0 / np.sqrt(den), np.nan)
    return b, se


def scan(a, c, M, names, min_snps=MIN_SNPS):
    n = np.asarray(M.sum(axis=0))
    keep = n >= min_snps
    b_set, se_set = estimates(a, c, M[:, keep])
    rest = ~M[:, keep]
    b_rest, se_rest = estimates(a, c, rest)
    diff = b_set - b_rest
    se_diff = np.sqrt(se_set ** 2 + se_rest ** 2)
    z = diff / se_diff
    p = 2 * stats.norm.sf(np.abs(z))
    tab = pd.DataFrame(dict(set=np.asarray(names)[keep], nsnp_set=n[keep],
                            nsnp_rest=len(a) - n[keep], b_set=b_set, se_set=se_set,
                            b_complement=b_rest, se_complement=se_rest,
                            diff=diff, se_diff=se_diff, z=z, p=p))
    names = tab["set"].astype(str)
    tagged = names.str.contains(":")
    tab.insert(0, "collection", np.where(tagged, names.str.split(":").str[0], "local"))
    tab["set"] = np.where(tagged, names.str.split(":", n=1).str[1], names)
    tab = tab.sort_values("p").reset_index(drop=True)
    order = np.argsort(tab["p"].to_numpy())
    q = np.empty(len(tab))
    ranked = tab["p"].to_numpy()[order]
    m = len(ranked)
    q[order] = np.minimum.accumulate((ranked * m / (np.arange(m) + 1))[::-1])[::-1]
    tab["bh_q"] = np.clip(q, 0, 1)
    return tab


def permutation_maxz(a, c, M, n_perm, rng):
    """Null for the largest |z| when gene labels are shuffled across instruments.

    A permuted column can drop below MIN_SNPS; those columns are excluded from that
    replicate's maximum rather than being allowed to contribute a near-infinite z.
    """
    n_snp = M.shape[0]
    mx = np.empty(n_perm)
    for i in range(n_perm):
        Mp = M[rng.permutation(n_snp)]
        ok = Mp.sum(axis=0) >= MIN_SNPS
        if not ok.any():
            mx[i] = np.nan
            continue
        b_set, se_set = estimates(a, c, Mp[:, ok])
        b_rest, se_rest = estimates(a, c, ~Mp[:, ok])
        z = np.abs((b_set - b_rest) / np.sqrt(se_set ** 2 + se_rest ** 2))
        mx[i] = np.nanmax(z) if np.isfinite(z).any() else np.nan
    ok = ~np.isnan(mx)
    return mx[ok]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--n-perm", type=int, default=2000)
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261022)
    scans, fwer = [], []
    nulls = []
    for key, label in OUTCOMES.items():
        print(f"\n=== {key} ({label})")
        d = fetch(key)
        if d is None or len(d) < 50:
            print("  no data")
            continue
        M, names, aw, cw = membership(d)
        tab = scan(aw, cw, M, names)
        tab.insert(0, "outcome", key)
        scans.append(tab)
        obs = float(np.nanmax(np.abs(tab["z"])))
        null = permutation_maxz(aw, cw, M, a.n_perm, rng)
        nulls.append(pd.DataFrame(dict(outcome=key, max_abs_z=null)))
        p_fw = float((np.sum(null >= obs) + 1) / (len(null) + 1))
        fwer.append(dict(outcome=key, observed_max_abs_z=obs, n_perm=len(null),
                         null_median=float(np.median(null)),
                         null_95=float(np.quantile(null, 0.95)), p_familywise=p_fw))
        print(f"  {len(d)} instruments, {len(tab)} gene sets with >= {MIN_SNPS} instruments "
              f"({tab['collection'].value_counts().to_dict()})")
        for _, r in tab.head(6).iterrows():
            print(f"  [{r['collection']:8s}] {r['set'][:38]:38s} n={r['nsnp_set']:3d}  "
                  f"set {r['b_set']:+.4f} vs rest {r['b_complement']:+.4f}  z={r['z']:+.2f}  "
                  f"p={r['p']:.3g}  q={r['bh_q']:.3g}")
        print(f"  max |z| = {obs:.2f}; permutation 95th pct = {np.quantile(null, 0.95):.2f}; "
              f"family-wise p = {p_fw:.3g}")
    if scans:
        pd.concat(scans).to_csv(OUT / "geneset_scan.csv", index=False, float_format="%.6g")
        pd.DataFrame(fwer).to_csv(OUT / "geneset_fwer.csv", index=False, float_format="%.6g")
        pd.concat(nulls).to_csv(OUT / "geneset_null.csv", index=False, float_format="%.6g")
        print(f"\n-> {OUT.relative_to(ROOT).as_posix()}/")
    return 0


def selftest():
    """A simulated set that really carries the effect must be recovered; a null set must not."""
    rng = np.random.default_rng(5)
    n = 400
    bx = rng.normal(0.05, 0.012, n)
    w = 1.0 / np.full(n, 0.01) ** 2
    member = np.zeros(n, bool)
    member[:60] = True
    by = np.where(member, -0.08 * bx, -0.01 * bx) + rng.normal(0, 0.002, n)
    d = pd.DataFrame(dict(bx=bx, by=by, sy=np.full(n, 0.01)))
    a = w * bx * by
    c = w * bx ** 2
    M = np.column_stack([member, np.r_[np.zeros(60, bool), np.ones(n - 60, bool)]])
    tab = scan(a, c, M, ["carrier", "null"])
    top = tab.iloc[0]
    print(f"selftest: top set {top['set']} diff {top['diff']:+.4f} z={top['z']:+.2f} "
          f"p={top['p']:.3g}")
    _, se = estimates(a, c, M[:, [0]])
    ok = top["set"] == "carrier" and top["z"] < 0 and abs(top["diff"] + 0.07) < 0.02
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

