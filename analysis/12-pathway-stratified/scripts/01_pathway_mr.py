"""Module 12 鈥?does the diffuse eBMD effect live in any bone pathway?

Module 09 showed the effect is diffuse at the *locus* level: 589 instruments map to 519
genes and only one locus (ABO) colocalises with heart failure. Diffuse loci do not imply
diffuse pathways 鈥?a hundred variants in WNT/osteoblast genes could still carry most of
the effect. This module asks that question directly and answers it in four ways:

  1. CONTRIBUTION DECOMPOSITION. For fixed-effect IVW the estimate is exactly
         beta = sum_i w_i bx_i by_i / sum_i w_i bx_i^2,   w_i = 1/sy_i^2
     so the effect splits additively across instrument groups: the share of group g is
         c_g = sum_{i in g} w_i bx_i by_i / sum_all w_i bx_i^2,   sum_g c_g = beta.
     This is an identity, not an approximation: it says where the estimate comes from.
  2. PER-GROUP ESTIMATE. IVW using only that group's instruments, with its own SE.
  3. BETWEEN-GROUP HETEROGENEITY. Cochran's Q across the group estimates.
  4. PERMUTATION NULL. Re-assign the group labels at random (same group sizes) and
     recompute Q 5,000 times. This is the test that pathway structure is real: if the
     curated partition produces no more heterogeneity than arbitrary partitions of the
     same size, the effect is diffuse in pathways as well as in loci.

Groups are pre-specified curated gene sets, applied in the listed order with first match
winning, so membership is reproducible and visible in the output tables. Anything not
listed is `unassigned` and is reported as its own group rather than dropped.

Outputs (results/tables/12-pathway-stratified/):
  pathway_contributions.csv     share of the effect and per-group IVW, per outcome
  pathway_group_estimates.csv   IVW within each group
  pathway_heterogeneity.csv     Q_between, its p, and the permutation p
  pathway_leave_one_out.csv     full-set estimate with each group removed

Usage:
  python analysis/12-pathway-stratified/scripts/01_pathway_mr.py
  python analysis/12-pathway-stratified/scripts/01_pathway_mr.py --self-test
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
EXPOSURE = "eBMD"
OUTCOMES = {"HF": "HERMES (overlapping)", "FG_HF": "FinnGen R12 (independent)"}

# Curated, mutually exclusive by first-match order. These are canonical bone-biology gene
# sets, not differential-expression output, so a gene's membership does not depend on
# these data and the partition cannot be tuned to the answer.
GROUPS = [
    ("osteoclast_resorption", [
        "TNFSF11", "TNFRSF11A", "TNFRSF11B", "CTSK", "ACP5", "CA2", "CLCN7", "OSTM1",
        "TCIRG1", "SNX10", "PLEKHM1", "SQSTM1", "CSF1", "CSF1R", "MITF", "NFATC1",
        "CALCR", "SRC", "ITGB3", "ATP6V0D2", "JDP2"]),
    ("wnt_osteoblast", [
        "SOST", "LRP4", "LRP5", "LRP6", "WNT16", "WNT4", "WNT5B", "WNT10B", "WLS",
        "DKK1", "DKK2", "SFRP1", "SFRP4", "AXIN1", "AXIN2", "GSK3B", "CTNNB1",
        "RUNX2", "SP7", "SOX6", "SOX9", "MSX2", "TWIST1", "TWIST2", "MEF2C",
        "PTHLH", "PTH1R", "GNAS", "GNA11", "FGFR1", "FGFRL1", "TGFB1", "BMP2", "BMP4",
        "BMP6", "BMP7", "BMPR1B", "BMPR2", "ACVR1", "ACVR2A", "SMAD3", "SMAD9",
        "NOG", "GREM1", "GREM2", "CHRD", "HDAC4", "HDAC5", "KDM6B", "EN1", "EN2",
        "GLI2", "IHH", "PTCH1", "SUFU", "ZBTB40", "CPED1", "FAM210A", "RSPO3",
        "DACT1", "PAPPA2", "CRB2"]),
    ("matrix_mineralisation", [
        "BGLAP", "MGP", "SPP1", "IBSP", "ALPL", "PHOSPHO1", "ENPP1", "ANKH", "SPARC",
        "PHEX", "DMP1", "COL1A1", "COL1A2", "COL2A1", "COL5A1", "COL6A1", "COL9A1",
        "COL11A1", "COL11A2", "COL27A1", "BGN", "DCN", "LUM", "SERPINH1", "CRTAP",
        "P3H1", "PLOD2", "MMP13", "VDR", "CYP24A1", "CYP27B1", "GC", "LRP2",
        "CUBN", "SLC34A1", "SLC34A3", "SLC9A3R1", "FGF23", "KL", "FAM20C"]),
    ("sex_hormone", [
        "ESR1", "ESR2", "AR", "SRD5A2", "CYP19A1", "SHBG", "LHCGR", "FSHR", "INHBA",
        "INHBB", "GREB1", "KLK3", "HSD17B2", "CYP11A1", "SRD5A1"]),
    ("adiposity_bodycomp", [
        "FTO", "MC4R", "BDNF", "SH2B1", "TMEM18", "SEC16B", "GNPDA2", "FAIM2",
        "KCTD15", "MTCH2", "NPC1", "POMC", "LEPR", "LEP", "PPARG", "ADRB3", "UCP1",
        "IRS1", "PPARGC1A", "ADCY3", "CREBBP", "PCSK1", "SIM1", "CADM1", "GIPR",
        "TFAP2B", "ETV5", "MAP2K5", "QPCTL", "NRXN3", "HOXB5", "SLC39A8"]),
    ("height_growth", [
        "FGFR3", "GDF5", "HMGA2", "ACAN", "NPR3", "IGF1", "IGF1R", "GH1", "GHR",
        "SOX5", "SHOX", "NPPB", "PRKG2", "ADAMTSL3"]),
    ("immune_inflammatory", [
        "IL6R", "IL6", "IL1B", "TNF", "CRP", "NLRP3", "IL17A", "TNFRSF1B",
        "HLA-DRB1", "HLA-B", "HLA-C", "HLA-DQA1", "ABO", "FUT2"]),
]


def group_of(gene):
    """First group (in the listed order) that contains the gene."""
    if not isinstance(gene, str) or not gene:
        return "unassigned"
    up = gene.upper()
    if up.startswith("HLA-"):
        return "immune_inflammatory"
    for name, genes in GROUPS:
        if up in genes:
            return name
    return "unassigned"


def fetch(outcome):
    """eBMD instruments aligned to the outcome: bx, sx, by, sy per SNP."""
    inst = mr.load_instruments(EXPOSURE)
    out = mr.outcome_at_snps(outcome, inst["SNP"].tolist())
    if out is None or out.empty:
        return None
    return mr.harmonise(inst, out)


def contributions(d, labels):
    """Exact fixed-effect split of the IVW estimate across labels."""
    bx = d["bx"].to_numpy(float)
    by = d["by"].to_numpy(float)
    w = 1.0 / d["sy"].to_numpy(float) ** 2
    denom = float(np.sum(w * bx ** 2))
    beta_fe = float(np.sum(w * bx * by)) / denom
    rows = []
    for g, idx in labels.groupby(labels).groups.items():
        m = d.index.isin(idx)
        num = float(np.sum(w[m] * bx[m] * by[m]))
        rows.append(dict(group=g, nsnp=int(m.sum()), contribution_b=num / denom,
                         share_pct=100 * num / denom / beta_fe if beta_fe else np.nan))
    return beta_fe, pd.DataFrame(rows)


def group_estimates(d, labels):
    """IVW within each group, for the between-group heterogeneity test."""
    rows = []
    for g, idx in labels.groupby(labels).groups.items():
        m = d.index.isin(idx)
        if m.sum() < 3:
            continue
        r = mr.ivw(d["bx"].to_numpy(float)[m], d["by"].to_numpy(float)[m],
                   d["sx"].to_numpy(float)[m], d["sy"].to_numpy(float)[m], random=True)
        rows.append(dict(group=g, nsnp=int(m.sum()), b=r["b"], se=r["se"], p=r["p"],
                         q_p=r["q_p"],
                         mean_F=float(np.mean((d["bx"].to_numpy(float)[m]
                                               / d["sx"].to_numpy(float)[m]) ** 2))))
    return pd.DataFrame(rows)


def between_q(est):
    if len(est) < 2:
        return np.nan, np.nan, np.nan
    w = 1.0 / est["se"].to_numpy(float) ** 2
    b = est["b"].to_numpy(float)
    pooled = float(np.sum(w * b) / np.sum(w))
    q = float(np.sum(w * (b - pooled) ** 2))
    return q, len(b) - 1, float(stats.chi2.sf(q, len(b) - 1))


def permute_labels(labels, sizes, rng):
    """Shuffle group labels at the SNP level, preserving the group sizes."""
    pool = np.concatenate([np.repeat(g, n) for g, n in sizes.items()])
    rng.shuffle(pool)
    return pd.Series(pool, index=labels.index)


def analyse(d, name, n_perm, rng):
    genes = pd.read_csv(GENES).set_index("SNP")["nearest_gene"]
    lab = pd.Series([group_of(genes.get(s, "")) for s in d["SNP"]], index=d.index, name="group")
    beta_fe, contrib = contributions(d, lab)
    est = group_estimates(d, lab)
    q, dfq, p_q = between_q(est)

    sizes = lab.value_counts().to_dict()
    null_q = np.empty(n_perm)
    for i in range(n_perm):
        e = group_estimates(d, permute_labels(lab, sizes, rng))
        null_q[i] = between_q(e)[0] if len(e) >= 2 else np.nan
    ok = ~np.isnan(null_q)
    p_perm = float((np.sum(null_q[ok] >= q) + 1) / (ok.sum() + 1)) if ok.any() else np.nan

    m = est.set_index("group")
    contrib["outcome"] = name
    contrib["group_b"] = contrib["group"].map(m["b"])
    contrib["group_se"] = contrib["group"].map(m["se"])
    contrib["group_p"] = contrib["group"].map(m["p"])
    contrib["mean_F"] = contrib["group"].map(m["mean_F"])
    contrib = contrib.sort_values("contribution_b", ascending=False).reset_index(drop=True)

    loo = []
    for g in contrib["group"]:
        keep = lab != g
        if keep.sum() < 10:
            continue
        r = mr.ivw(d["bx"].to_numpy(float)[keep], d["by"].to_numpy(float)[keep],
                   d["sx"].to_numpy(float)[keep], d["sy"].to_numpy(float)[keep], random=True)
        loo.append(dict(outcome=name, removed_group=g, nsnp=int(keep.sum()), b=r["b"],
                        se=r["se"], p=r["p"], delta_b=r["b"] - beta_fe,
                        pct_change=100 * (r["b"] - beta_fe) / beta_fe if beta_fe else np.nan))
    het = dict(outcome=name, beta_fixed=beta_fe, n_groups=int(len(est)),
               q_between=q, df=dfq, p_q=p_q, p_permutation=p_perm,
               n_perm=int(ok.sum()), null_q_median=float(np.nanmedian(null_q[ok])))
    return contrib, pd.DataFrame(loo), het, est.assign(outcome=name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--n-perm", type=int, default=5000)
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261022)
    cs, ls, hs, es = [], [], [], []
    for key, label in OUTCOMES.items():
        print(f"\n=== {key} ({label})")
        d = fetch(key)
        if d is None or len(d) < 50:
            print("  no data")
            continue
        contrib, loo, het, est = analyse(d, key, a.n_perm, rng)
        cs.append(contrib)
        ls.append(loo)
        hs.append(het)
        es.append(est)
        print(f"  {len(d)} instruments, fixed-effect beta = {het['beta_fixed']:+.4f}")
        for _, r in contrib.iterrows():
            print(f"  {r['group']:22s} n={r['nsnp']:4d}  contribution "
                  f"{r['contribution_b']:+.4f} ({r['share_pct']:+6.1f}%)  "
                  f"group IVW {r['group_b']:+.4f} (p={r['group_p']:.3g})")
        print(f"  between-group Q = {het['q_between']:.2f} (df {het['df']}, p={het['p_q']:.3g}); "
              f"permutation p = {het['p_permutation']:.3g} over {het['n_perm']} shuffles")
    if cs:
        pd.concat(cs).to_csv(OUT / "pathway_contributions.csv", index=False, float_format="%.6g")
        pd.concat(ls).to_csv(OUT / "pathway_leave_one_out.csv", index=False, float_format="%.6g")
        pd.DataFrame(hs).to_csv(OUT / "pathway_heterogeneity.csv", index=False,
                                float_format="%.6g")
        pd.concat(es).to_csv(OUT / "pathway_group_estimates.csv", index=False,
                             float_format="%.6g")
        print(f"\n-> {OUT.relative_to(ROOT).as_posix()}/")
    return 0


def selftest():
    """Two simulated groups: one carries the effect, the other is pure noise.

    The decomposition must sum to the fixed-effect estimate, and the between-group test
    must separate a real split from noise.
    """
    rng = np.random.default_rng(11)
    n1 = n2 = 300
    bx = rng.normal(0.05, 0.012, n1 + n2)
    sx = np.full(n1 + n2, 0.005)
    sy = np.full(n1 + n2, 0.01)
    by = np.concatenate([-0.05 * bx[:n1] + rng.normal(0, 0.0005, n1),
                         rng.normal(0, 0.002, n2)])
    d = pd.DataFrame(dict(bx=bx, by=by, sx=sx, sy=sy))
    lab = pd.Series(["carrier"] * n1 + ["null"] * n2, index=d.index)
    beta_fe, contrib = contributions(d, lab)
    if abs(contrib["contribution_b"].sum() - beta_fe) > 1e-9:
        print(f"FAIL: shares sum to {contrib['contribution_b'].sum():.6f}, beta {beta_fe:.6f}")
        return 1
    est = group_estimates(d, lab)
    q, _, p_q = between_q(est)
    share = contrib.loc[contrib.group == "carrier", "share_pct"].iloc[0]
    print(f"selftest: beta_FE {beta_fe:+.4f}, carrier group carries {share:.1f}% of it")
    print(f"selftest: between-group Q = {q:.2f} (p={p_q:.3g})")
    ok = share > 55 and p_q < 0.01
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

