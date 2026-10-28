"""Module 07-scrna 鈥?where in the calcified plaque does the bone programme live?

GSE159677 is single-cell RNA-seq of carotid endarterectomy tissue from three patients,
sequenced as whole calcified atherosclerotic core (AC) and patient-matched proximal
adjacent (PA) artery. It is the only dataset in the project that can locate a cell type
rather than a bulk average.

Three questions:
  1. Which cell types are present, and how does their composition differ between the
     calcified core and the adjacent tissue?
  2. Which cell type carries an osteogenic / osteoblast-like programme?
  3. Which cell type expresses the bone-pathway genes we identified, i.e. where in the
     plaque does the BMD signal actually land?

Library assignment. The aggregate matrix shipped by GEO (`AGGREGATEMAPPED`) carries
barcodes suffixed -1..-6 but no sample label, and the per-sample bundles inside RAW.tar
cannot be downloaded individually (they 404). The assignment was therefore *verified*
rather than assumed: the one per-sample bundle obtained (GSM4837523, patient 1 core)
contains 11,015 barcodes, and exactly 11,015 aggregate barcodes carry the -2 suffix with a
100 % sequence match, while every other suffix matches 10-48 of them by chance. Cell counts
then follow library size across all six (3.4k-16k cells, 27-92 MB), and the core libraries
are consistently the larger ones. Core = {-2, -4, -6}, adjacent = {-1, -3, -5}.

Outputs (results/tables/07-scrna/):
  celltype_composition.csv   cell-type fractions in core vs adjacent
  programme_scores.csv       per-cell-type programme means
  osteo_by_tissue.csv        osteogenic/calcification scores, core vs adjacent
  markers_used.csv           marker sets that were actually measurable

Usage:
  python analysis/07-scrna/scripts/01_scrna_localise.py
  python analysis/07-scrna/scripts/01_scrna_localise.py --self-test
"""
import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import TruncatedSVD

ROOT = Path(__file__).resolve().parents[3]
MAT = ROOT / "data" / "processed" / "scrna" / "GSE159677" / "outs" / "filtered_feature_bc_matrix"
OUT = ROOT / "results" / "tables" / "07-scrna"

CORE_SUFFIX = {"2", "4", "6"}
ADJ_SUFFIX = {"1", "3", "5"}

MARKERS = {
    "smooth muscle": ["ACTA2", "TAGLN", "MYH11", "CNN1", "DES"],
    "fibroblast": ["LUM", "DCN", "PDGFRA", "COL1A2", "COL3A1"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "EGFL7"],
    "macrophage": ["CD68", "CD14", "LYZ", "MRC1", "C1QA"],
    "T cell": ["CD3D", "CD3E", "IL7R", "CD2"],
    "B cell / plasma": ["MS4A1", "CD79A", "IGHG1", "MZB1"],
    "mast": ["TPSAB1", "CPA3", "KIT"],
    "osteogenic": ["RUNX2", "SP7", "ALPL", "BGLAP", "SPP1", "MMP13", "COL1A1", "SOX9", "IBSP"],
    "chondrogenic": ["SOX9", "ACAN", "COL2A1", "COMP"],
    "calcification": ["ALPL", "ENPP1", "MGP", "SPP1", "BMP2"],
    "inflammation": ["IL6", "IL1B", "TNF", "CCL2", "NLRP3"],
}

CELLTYPE_MARKERS = ["smooth muscle", "fibroblast", "endothelial", "macrophage",
                    "T cell", "B cell / plasma", "mast"]

# bone-pathway genes: the loci the BMD instruments sit in plus the mediators tested in
# module 06, used as a fixed, pre-specified list rather than selected from the data
BMD_GENES = ["ESR1", "WNT16", "CPED1", "LRP5", "SOST", "RUNX2", "SP7", "ALPL",
             "BGLAP", "COL1A1", "VDR", "PTH1R", "FGF23", "DKK1", "MGP"]


def load_matrix():
    """Return (cells x genes CSR float32, gene symbols, barcodes)."""
    from scipy.io import mmread
    with gzip.open(MAT / "matrix.mtx.gz", "rb") as f:
        m = mmread(f).tocsr()                       # genes x cells
    with gzip.open(MAT / "features.tsv.gz", "rt") as f:
        genes = [line.split("\t")[1] for line in f]
    with gzip.open(MAT / "barcodes.tsv.gz", "rt") as f:
        barcodes = [line.strip() for line in f]
    return m.T.tocsr().astype(np.float32), genes, barcodes


def qc_and_normalise(X, min_umi=200, min_cells=10):
    """Drop low-UMI cells and rarely detected genes, then CP10k + log1p."""
    lib = np.asarray(X.sum(axis=1)).ravel()
    gene_cells = np.asarray((X > 0).sum(axis=0)).ravel()
    keep_cells = lib >= min_umi
    keep_genes = gene_cells >= min_cells
    X = X[keep_cells][:, keep_genes]
    lib = np.asarray(X.sum(axis=1)).ravel()
    X = sp.diags(1e4 / lib) @ X
    X.data = np.log1p(X.data)
    return X.tocsr(), keep_cells, keep_genes


def score_sets(X, genes, sets):
    """Mean log-expression of each gene set, over the genes actually measurable."""
    idx = {g: i for i, g in enumerate(genes)}
    out, used = {}, {}
    for name, gl in sets.items():
        cols = [idx[g] for g in gl if g in idx]
        used[name] = [genes[c] for c in cols]
        out[name] = (np.asarray(X[:, cols].mean(axis=1)).ravel().astype(np.float32)
                     if cols else np.zeros(X.shape[0], np.float32))
    return out, used


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--clusters", type=int, default=14)
    ap.add_argument("--hvg", type=int, default=2000)
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    print("loading matrix ...")
    X, genes, barcodes = load_matrix()
    print(f"  raw {X.shape[0]:,} cells x {X.shape[1]:,} genes, {X.nnz:,} non-zero")
    X, okc, okg = qc_and_normalise(X)
    genes = [g for g, k in zip(genes, okg) if k]
    barcodes = [b for b, k in zip(barcodes, okc) if k]
    print(f"  after QC {X.shape[0]:,} cells x {X.shape[1]:,} genes")

    suffix = np.array([b.rsplit("-", 1)[-1] for b in barcodes])
    tissue = np.where(np.isin(suffix, list(CORE_SUFFIX)), "core",
                      np.where(np.isin(suffix, list(ADJ_SUFFIX)), "adjacent", "unknown"))
    print(f"  tissue: core {int((tissue == 'core').sum()):,}, "
          f"adjacent {int((tissue == 'adjacent').sum()):,}, unknown {int((tissue == 'unknown').sum()):,}")
    keep = tissue != "unknown"
    X = X[keep]
    genes_arr = np.array(genes)
    barcodes = [b for b, k in zip(barcodes, keep) if k]
    suffix, tissue = suffix[keep], tissue[keep]

    print("feature selection + embedding ...")
    mean = np.asarray(X.mean(axis=0)).ravel()
    var = np.asarray(X.multiply(X).mean(axis=0)).ravel() - mean ** 2
    disp = np.where(var > 0, var / (mean + 1e-9), 0)
    hvg = np.argsort(disp)[::-1][: a.hvg]
    svd = TruncatedSVD(n_components=30, random_state=0)
    emb = svd.fit_transform(X[:, hvg])
    print(f"  {a.hvg} HVGs -> 30 components, {svd.explained_variance_ratio_.sum():.1%} variance")

    km = MiniBatchKMeans(n_clusters=a.clusters, random_state=0, batch_size=4096, n_init=5)
    df = pd.DataFrame({"cluster": km.fit_predict(emb)})

    print("scoring marker programmes ...")
    sc, used = score_sets(X, list(genes_arr), MARKERS)
    pd.DataFrame([(k, ";".join(v)) for k, v in used.items()],
                 columns=["programme", "genes_used"]).to_csv(OUT / "markers_used.csv", index=False)
    for k, v in sc.items():
        df[k] = v
    df["tissue"] = tissue
    df["suffix"] = suffix

    prof = df.groupby("cluster")[CELLTYPE_MARKERS].mean()
    overall = prof.values.mean()
    names = {c: (prof.loc[c].idxmax() if prof.loc[c].max() > overall else "unassigned")
             for c in prof.index}
    df["celltype"] = df["cluster"].map(names)
    print("  cluster -> cell type:", {int(k): v for k, v in sorted(names.items())})

    # ---- 1. composition
    comp = df.groupby(["celltype", "tissue"]).size().unstack(fill_value=0)
    comp_pct = (comp.T / comp.sum(axis=1)).T * 100
    comp_pct["n_cells"] = df["celltype"].value_counts()
    comp_pct.round(3).to_csv(OUT / "celltype_composition.csv")
    print("\ncell-type composition (% of cells per tissue):")
    print(comp_pct.round(2).to_string())

    # ---- 2. programme scores, core vs adjacent, per cell type
    rows = []
    for ct in sorted(df["celltype"].unique()):
        sub = df[df["celltype"] == ct]
        for prog in ["osteogenic", "chondrogenic", "calcification", "inflammation"]:
            c = sub.loc[sub.tissue == "core", prog]
            p = sub.loc[sub.tissue == "adjacent", prog]
            if len(c) < 20 or len(p) < 20:
                continue
            u = stats.mannwhitneyu(c, p, alternative="two-sided")
            rows.append(dict(celltype=ct, programme=prog, n_core=len(c), n_adjacent=len(p),
                             mean_core=float(c.mean()), mean_adjacent=float(p.mean()),
                             diff=float(c.mean() - p.mean()), p=float(u.pvalue)))
    prog_df = pd.DataFrame(rows)
    if len(prog_df):
        prog_df["p_adj"] = bh(prog_df["p"])
        prog_df.sort_values("p").round(6).to_csv(OUT / "osteo_by_tissue.csv", index=False)
        print("\nprogrammes, core vs adjacent (Mann-Whitney):")
        print(prog_df.sort_values(["programme", "p"])
              [["celltype", "programme", "mean_core", "mean_adjacent", "diff", "p", "p_adj"]]
              .round(4).to_string(index=False))

    # ---- 3. bone-pathway genes per cell type
    idx = {g: i for i, g in enumerate(genes_arr)}
    present = [g for g in BMD_GENES if g in idx]
    print(f"\nbone-pathway genes measurable: {len(present)}/{len(BMD_GENES)}: {', '.join(present)}")
    df["bmd_gene_score"] = (np.asarray(X[:, [idx[g] for g in present]].mean(axis=1)).ravel()
                            if present else 0.0)
    per_ct = df.groupby("celltype").agg(
        n=("bmd_gene_score", "size"),
        bmd_gene_score=("bmd_gene_score", "mean"),
        osteogenic=("osteogenic", "mean"),
        calcification=("calcification", "mean"),
        inflammation=("inflammation", "mean")).sort_values("bmd_gene_score", ascending=False)
    per_ct.round(4).to_csv(OUT / "programme_scores.csv")
    print("\nper cell type (sorted by bone-pathway gene score):")
    print(per_ct.round(3).to_string())

    # per-cell labels, so module 09 can join instrument genes to cell types without
    # repeating the clustering
    meta = pd.DataFrame({"barcode": barcodes, "celltype": df["celltype"].to_numpy(),
                         "tissue": df["tissue"].to_numpy(), "suffix": df["suffix"].to_numpy(),
                         "cluster": df["cluster"].to_numpy()})
    meta.to_csv(OUT / "cell_metadata.csv.gz", index=False, compression="gzip")
    print(f"\n-> {OUT / 'cell_metadata.csv.gz'} ({len(meta):,} cells)")
    return 0


def selftest():
    """Two planted populations; the scoring must recover the osteogenic programme in the
    one that carries it and remain flat in the other."""
    rng = np.random.default_rng(7)
    genes = ["ACTA2", "MYH11", "CD68", "LYZ", "RUNX2", "ALPL", "SPP1", "BGLAP",
             "PECAM1", "VWF", "COL1A1", "SOX9"] + [f"G{i}" for i in range(300)]
    n1, n2 = 400, 400
    X = np.zeros((n1 + n2, len(genes)), np.float32)
    X[:n1, [0, 1, 10]] = rng.gamma(2, 1.0, (n1, 3))
    X[n1:, [2, 3]] = rng.gamma(2, 1.0, (n2, 2))
    X[n1:, [4, 5, 6, 7, 11]] = rng.gamma(3, 1.2, (n2, 5))
    X[:, 8:10] = rng.gamma(1, 0.3, (n1 + n2, 2))
    X[:, 12:] = rng.gamma(1, 0.2, (n1 + n2, 300))
    Xn, _, okg = qc_and_normalise(sp.csr_matrix(X), min_umi=10, min_cells=5)
    sc, _ = score_sets(Xn, [g for g, k in zip(genes, okg) if k], MARKERS)
    osteo_hi, osteo_lo = sc["osteogenic"][n1:].mean(), sc["osteogenic"][:n1].mean()
    sm_hi, sm_lo = sc["smooth muscle"][:n1].mean(), sc["smooth muscle"][n1:].mean()
    print(f"selftest: osteogenic {osteo_hi:.3f} (planted) vs {osteo_lo:.3f}")
    print(f"selftest: smooth muscle {sm_hi:.3f} vs {sm_lo:.3f}")
    ok = osteo_hi > osteo_lo * 1.5 and sm_hi > sm_lo * 1.5
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

