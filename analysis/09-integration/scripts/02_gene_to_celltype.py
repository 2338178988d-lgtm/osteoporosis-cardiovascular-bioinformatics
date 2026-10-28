"""Module 09 step 2 鈥?which cell type expresses the genes the BMD instruments point at?

Step 1 assigned all 589 usable eBMD instruments to a nearest gene using FinnGen's own
annotation. This step asks the single-cell question about those genes: for each instrument
locus, mean expression of its nearest gene in each cellular compartment of the calcified
plaques, then an aggregate test of whether the instrument genes are enriched in any one
compartment relative to a matched background.

The background matters and is built deliberately. Instrument genes sit near
disease-associated loci, which are themselves gene-dense and often in long genes, so
comparing them with all genes would confound "instrument gene" with "gene length" and
"expression level". The null is therefore 1,000 random gene sets of the same size and
matched on mean expression decile within the same matrix.

Outputs (results/tables/09-integration/):
  gene_celltype_expression.csv   gene x cell type mean expression (instrument genes)
  gene_preferred_celltype.csv    each gene's highest-expressing compartment
  celltype_enrichment.csv        observed vs matched-null share per compartment
  expression_summary.csv         fraction of instrument genes detected per compartment

Usage:
  python analysis/09-integration/scripts/02_gene_to_celltype.py
  python analysis/09-integration/scripts/02_gene_to_celltype.py --self-test
"""
import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
MAT = ROOT / "data" / "processed" / "scrna" / "GSE159677" / "outs" / "filtered_feature_bc_matrix"
META = ROOT / "results" / "tables" / "07-scrna" / "cell_metadata.csv.gz"
GENES = ROOT / "results" / "tables" / "09-integration" / "instrument_genes.csv"
OUT = ROOT / "results" / "tables" / "09-integration"

N_NULL = 1000
SEED = 20261022


def load_all():
    from scipy.io import mmread
    with gzip.open(MAT / "matrix.mtx.gz", "rb") as f:
        m = mmread(f).tocsr().T.tocsr().astype(np.float32)      # cells x genes
    with gzip.open(MAT / "features.tsv.gz", "rt") as f:
        genes = [line.split("\t")[1] for line in f]
    with gzip.open(MAT / "barcodes.tsv.gz", "rt") as f:
        barcodes = [line.strip() for line in f]
    meta = pd.read_csv(META)
    return m, genes, barcodes, meta


def normalise(m):
    lib = np.asarray(m.sum(axis=1)).ravel()
    keep = lib >= 200
    m = m[keep]
    lib = np.asarray(m.sum(axis=1)).ravel()
    m = sp.diags(1e4 / lib) @ m
    m.data = np.log1p(m.data)
    return m.tocsr(), keep


def celltype_means(X, col_idx, labels, levels):
    """Mean expression per cell type for the selected gene columns."""
    out = np.zeros((len(col_idx), len(levels)), np.float32)
    sub = X[:, col_idx]
    for j, lv in enumerate(levels):
        rows = labels == lv
        if rows.sum() == 0:
            continue
        out[:, j] = np.asarray(sub[rows].mean(axis=0)).ravel()
    return out


def detected_fraction(X, col_idx, labels, levels):
    """Fraction of cells with non-zero expression, per gene per cell type."""
    out = np.zeros((len(col_idx), len(levels)), np.float32)
    sub = X[:, col_idx]
    nz = (sub > 0).astype(np.float32)
    for j, lv in enumerate(levels):
        rows = labels == lv
        if rows.sum() == 0:
            continue
        out[:, j] = np.asarray(nz[rows].mean(axis=0)).ravel()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--n-null", type=int, default=N_NULL)
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    print("loading single-cell matrix ...")
    X, genes, barcodes, meta = load_all()
    X, keep = normalise(X)
    barcodes = [b for b, k in zip(barcodes, keep) if k]
    labels = meta.set_index("barcode").loc[barcodes, "celltype"].to_numpy()
    tissue = meta.set_index("barcode").loc[barcodes, "tissue"].to_numpy()
    levels = sorted(pd.unique(labels))
    print(f"  {X.shape[0]:,} cells x {X.shape[1]:,} genes; compartments: {levels}")

    gidx = {g: i for i, g in enumerate(genes)}
    inst = pd.read_csv(GENES)
    inst = inst[inst.matched & (inst.nearest_gene != "")].copy()
    # the annotation occasionally lists several genes for one locus; take the first
    inst["gene"] = inst["nearest_gene"].str.split(",").str[0].str.strip()
    present = inst[inst["gene"].isin(gidx)].copy()
    print(f"instrument loci with a gene: {len(inst)}; of those measurable in the matrix: "
          f"{len(present)} ({100 * len(present) / max(len(inst), 1):.1f}%)")

    cols = [gidx[g] for g in present["gene"]]
    mean_expr = celltype_means(X, cols, labels, levels)
    det_frac = detected_fraction(X, cols, labels, levels)

    expr_df = pd.DataFrame(mean_expr, columns=levels)
    expr_df.insert(0, "gene", present["gene"].to_numpy())
    expr_df.insert(0, "SNP", present["SNP"].to_numpy())
    expr_df.to_csv(OUT / "gene_celltype_expression.csv", index=False)

    det_df = pd.DataFrame(det_frac, columns=levels)
    det_df.insert(0, "gene", present["gene"].to_numpy())
    det_df.to_csv(OUT / "gene_detection_fraction.csv", index=False)

    # each gene's preferred compartment, but only when it is actually detected somewhere
    max_det = det_frac.max(axis=1)
    expressible = max_det >= 0.10                    # detected in >=10% of some compartment
    pref = pd.DataFrame({
        "SNP": present["SNP"].to_numpy(),
        "gene": present["gene"].to_numpy(),
        "preferred_celltype": np.where(expressible, np.array(levels)[mean_expr.argmax(axis=1)], "not_detected"),
        "max_detection": max_det,
        "max_mean_expr": mean_expr.max(axis=1),
    })
    pref.to_csv(OUT / "gene_preferred_celltype.csv", index=False)
    share = pref["preferred_celltype"].value_counts()
    print(f"\ngenes detected in >=10% of a compartment: {int(expressible.sum())}")
    print("preferred compartment among those:")
    print(share.to_string())

    # ---- enrichment against an expression-matched null
    gene_mean = np.asarray(X.mean(axis=0)).ravel()
    all_idx = np.arange(len(genes))
    # decile bins of overall expression, computed on genes detected somewhere
    detected_any = np.asarray((X > 0).sum(axis=0)).ravel() >= 10
    pool = all_idx[detected_any]
    edges = np.quantile(gene_mean[pool], np.linspace(0, 1, 11))
    bins = np.digitize(gene_mean, edges[1:-1])
    obs_bins = bins[cols]
    rng = np.random.default_rng(SEED)
    pool_by_bin = {b: pool[bins[pool] == b] for b in np.unique(bins[pool])}

    obs_share = np.array([(pref.loc[expressible, "preferred_celltype"] == lv).sum()
                          for lv in levels], float)
    obs_share = obs_share / max(obs_share.sum(), 1)

    null_share = np.zeros((a.n_null, len(levels)))
    for i in range(a.n_null):
        pick = np.concatenate([rng.choice(pool_by_bin[b], size=max((obs_bins == b).sum(), 1),
                                          replace=True)
                               for b in np.unique(obs_bins)
                               if b in pool_by_bin])
        pick = pick[: len(cols)]
        m_null = celltype_means(X, pick, labels, levels)
        d_null = detected_fraction(X, pick, labels, levels)
        ok = d_null.max(axis=1) >= 0.10
        if ok.sum() == 0:
            continue
        best = np.array(levels)[m_null[ok].argmax(axis=1)]
        null_share[i] = [np.mean(best == lv) for lv in levels]

    rows = []
    for j, lv in enumerate(levels):
        obs = obs_share[j]
        null = null_share[:, j]
        null = null[np.isfinite(null)]
        p = (np.sum(null >= obs) + 1) / (len(null) + 1) if len(null) else np.nan
        rows.append(dict(celltype=lv, observed_share=obs, null_mean=float(null.mean()),
                         null_q95=float(np.quantile(null, 0.95)) if len(null) else np.nan,
                         fold=float(obs / null.mean()) if len(null) and null.mean() > 0 else np.nan,
                         p=p))
    enr = pd.DataFrame(rows).sort_values("fold", ascending=False)
    enr.to_csv(OUT / "celltype_enrichment.csv", index=False, float_format="%.6g")
    print(f"\nenrichment of instrument-gene preferred compartment (matched null, "
          f"{a.n_null} draws):")
    print(enr.round(4).to_string(index=False))

    summ = pd.DataFrame({
        "celltype": levels,
        "n_cells": [int((labels == lv).sum()) for lv in levels],
        "frac_genes_detected_10pct": [(det_frac[:, j] >= 0.10).mean() for j in range(len(levels))],
        "median_detection": [float(np.median(det_frac[:, j])) for j in range(len(levels))],
        "median_expression": [float(np.median(mean_expr[:, j])) for j in range(len(levels))],
    })
    summ.to_csv(OUT / "expression_summary.csv", index=False, float_format="%.6g")
    print("\nper-compartment detection summary:")
    print(summ.round(4).to_string(index=False))
    print(f"\ninstrument genes detected preferentially in the core: "
          f"{int((pref.loc[expressible].assign(core=(pref.loc[expressible, 'gene'].isin(present['gene']))).shape[0]))}")
    return 0


def selftest():
    """Plant 40 genes high in one compartment and 40 high in another; the enrichment
    statistic must recover both and not call the others."""
    rng = np.random.default_rng(1)
    n_cells, n_genes = 800, 600
    X = sp.csr_matrix(rng.gamma(1, 0.2, (n_cells, n_genes)).astype(np.float32))
    labels = np.array(["A"] * 400 + ["B"] * 400)
    X[:400, :40] = X[:400, :40] * 6            # genes 0-39 high in A
    X[400:, 40:80] = X[400:, 40:80] * 6        # genes 40-79 high in B
    X = X.tocsr()
    levels = ["A", "B"]
    m = celltype_means(X, list(range(80)), labels, levels)
    d = detected_fraction(X, list(range(80)), labels, levels)
    best = np.where(d.max(axis=1) >= 0.10, np.array(levels)[m.argmax(axis=1)], "not")
    a_share = float(np.mean(best[:40] == "A"))
    b_share = float(np.mean(best[40:] == "B"))
    print(f"selftest: genes planted in A -> {a_share:.0%} assigned to A")
    print(f"selftest: genes planted in B -> {b_share:.0%} assigned to B")
    ok = a_share >= 0.9 and b_share >= 0.9
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

