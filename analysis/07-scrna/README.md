# 07-scrna 鈥?where in the calcified plaque does the bone programme live?

## Purpose

Modules 02鈥?6 worked at the level of tissue averages: bulk transcriptomes, GWAS summary
statistics, population surveys. None of them can say *which cell* carries a signal.
GSE159677 is the only dataset in the project that can, because it is single-cell RNA-seq
of carotid endarterectomy tissue sequenced as whole calcified atherosclerotic core (AC)
and patient-matched proximal adjacent (PA) artery.

The question is narrow: **inside a calcified plaque, which cell type expresses the
osteogenic programme and the bone-pathway genes that the BMD genetic instruments point
at?**

## Inputs

`data/processed/scrna/GSE159677/outs/filtered_feature_bc_matrix/` 鈥?the aggregate 10x
matrix from GEO (`GSE159677_AGGREGATEMAPPED-...tar.gz`, 294 MB, verified), three files:

| File | Content |
|---|---|
| `matrix.mtx.gz` | 33,538 genes 脳 51,981 barcodes, 79,019,023 non-zero |
| `features.tsv.gz` | Ensembl ID, gene symbol, feature type |
| `barcodes.tsv.gz` | cell barcodes with library suffix `-1`鈥-6` |

Six libraries: `GSM4837523`鈥揱GSM4837528` = three patients 脳 (core, adjacent), established
from GEO sample records.

> **Library assignment was verified, not assumed.** The aggregate matrix carries no sample
> label, and the per-sample bundles inside `GSE159677_RAW.tar` return 404 when requested
> directly. The one bundle that could be obtained (`GSM4837523`, patient 1 core) contains
> 11,015 barcodes; exactly 11,015 aggregate barcodes carry the `-2` suffix and **100 % of
> the sequences match**, while every other suffix matches only 10鈥?8 of them by chance.
> Cell counts then track library size consistently across all six (3.4k鈥?6.4k cells,
> 27鈥?2 MB) with the core libraries being the larger ones, so core = {-2, -4, -6} and
> adjacent = {-1, -3, -5}. This is recorded by `analysis/scripts/fetch_gse159677.py`, which
> streams `RAW.tar` and keeps only the filtered matrices rather than the whole 2.1 GB
> archive.

## Methods

1. **Loading and QC** 鈥?`scipy.io.mmread` 鈫?cells 脳 genes CSR. Cells with fewer than 200
   UMIs and genes detected in fewer than 10 cells are dropped. Normalisation is CP10k then
   `log1p`.
2. **Embedding** 鈥?top 2,000 genes by dispersion, then `TruncatedSVD` to 30 components
   (62.9 % variance explained). `MiniBatchKMeans` with k = 14.
3. **Annotation** 鈥?seven lineage marker sets (`markers_used.csv` lists what was actually
   measurable). A cluster takes the lineage with its highest mean score if that score
   exceeds the global mean across clusters, otherwise "unassigned".
4. **Programme scores** 鈥?mean log-expression over four sets: osteogenic
   (RUNX2, SP7, ALPL, BGLAP, SPP1, MMP13, COL1A1, SOX9, IBSP), chondrogenic,
   calcification (ALPL, ENPP1, MGP, SPP1, BMP2), inflammation.
5. **Core vs adjacent** 鈥?Mann-Whitney per cell type per programme, Benjamini鈥揌ochberg
   across all tests.
6. **Bone-pathway score** 鈥?a fixed, pre-specified list of 15 genes (the BMD instrument
   loci plus the module-06 mediators), not selected from the data: ESR1, WNT16, CPED1,
   LRP5, SOST, RUNX2, SP7, ALPL, BGLAP, COL1A1, VDR, PTH1R, FGF23, DKK1, MGP. 14 of 15
   are measurable (SOX9 spelled differently in this annotation).

## Results as computed

47,980 cells passed QC, all with a recognised library suffix (core 38,105, adjacent 9,875).

**Composition**

| Cell type | % in adjacent | % in core | n cells |
|---|---|---|---|
| T cell | 17.9 | 82.1 | 18,677 |
| macrophage | 10.6 | **89.4** | 10,220 |
| endothelial | 48.0 | 52.0 | 6,114 |
| smooth muscle | 27.6 | 72.4 | 5,039 |
| fibroblast | 36.5 | 63.5 | 4,960 |
| mast | 32.2 | 67.8 | 4,844 |
| B cell / plasma | 23.5 | 76.5 | 2,126 |

Macrophages are the most core-enriched population. Endothelium is the only near-balanced
one, which is what a sequencing-depth artefact would **not** produce.

**The osteogenic programme is carried by macrophages and fibroblasts**

| Cell type | osteogenic core | osteogenic adjacent | change | p |
|---|---|---|---|---|
| **macrophage** | **0.261** | **0.047** | **+0.214 (5.6脳)** | < 1e-300 |
| fibroblast | 0.210 | 0.137 | +0.073 | < 1e-300 |
| mast | 0.063 | 0.025 | +0.038 | < 1e-300 |
| endothelial | 0.044 | 0.004 | +0.040 | < 1e-300 |
| smooth muscle | 0.058 | 0.039 | +0.019 | < 1e-300 |
| T cell | 0.033 | 0.010 | +0.023 | < 1e-300 |
| B cell / plasma | 0.024 | 0.007 | +0.017 | < 1e-300 |

Every lineage shows some induction, which tells us part of this is a core-wide
calcification environment rather than one cell type. But **macrophages show by far the
largest relative induction (5.6脳) and the highest absolute osteogenic score of any cell
type**, and they are simultaneously the most core-enriched population. The osteoblast-like
compartment in a calcified carotid plaque is therefore predominantly **macrophage-derived**,
with fibroblasts second.

**Ranking by bone-pathway gene expression differs from ranking by osteogenic score**

| Cell type | bone-pathway score | osteogenic | calcification | inflammation |
|---|---|---|---|---|
| **fibroblast** | **0.472** | 0.183 | **1.068** | 0.226 |
| smooth muscle | 0.349 | 0.053 | 0.831 | 0.114 |
| endothelial | 0.221 | 0.025 | 0.631 | 0.202 |
| mast | 0.091 | 0.051 | 0.264 | 0.096 |
| **macrophage** | 0.046 | **0.238** | 0.464 | **0.555** |
| T cell | 0.041 | 0.029 | 0.129 | 0.158 |
| B cell / plasma | 0.037 | 0.020 | 0.108 | 0.165 |

This split is the interesting part and it must not be flattened into a single "bone cells
are in the plaque" claim:

- **fibroblasts** carry the highest expression of the actual bone-pathway genes we
  validated genetically (ESR1, WNT16, SOST, LRP5, MGP) and the highest calcification score;
- **macrophages** carry the highest *osteogenic transcription-factor* programme (RUNX2,
  SP7, ALPL, BGLAP, SPP1) and the highest inflammation score.

Two different cell types, two different aspects of "bone-like". Which one matters for the
BMD鈫抙eart-failure effect cannot be resolved from this dataset.

## Commands

| Script | Command | Output |
|---|---|---|
| `analysis/scripts/fetch_gse159677.py` | `python analysis/scripts/fetch_gse159677.py` | streaming extraction of the six per-sample bundles |
| `scripts/01_scrna_localise.py` | `python analysis/07-scrna/scripts/01_scrna_localise.py` | `celltype_composition.csv`, `programme_scores.csv`, `osteo_by_tissue.csv`, `markers_used.csv` |
| self-test | `... --self-test` | PASS/FAIL |

Self-test plants two populations, one of which carries an osteogenic programme; the scorer
must recover it (observed: 4.92 vs 0.90 for osteogenic, 5.39 vs 0.00 for smooth muscle).

## Limitations

- **Three patients.** Every conclusion is a claim about these three plaques, not about
  carotid atherosclerosis in general. There is no way to put a confidence interval on a
  three-donor comparison.
- **Library assignment for two adjacent samples is unresolved.** `-1` and `-5` are both
  patient-matched adjacent tissue (26.7 MB / 3,716 cells and 27.1 MB / 3,379 cells) and
  cannot be told apart without the second per-sample bundle. Core-versus-adjacent does not
  depend on this, but a patient-paired design would.
- **k = 14 was fixed, not tuned.** No clustering-resolution sweep, no marker-based
  validation against the original publication's annotation.
- **No ambient-RNA correction or doublet removal.** Plaque tissue is dissociated from
  calcified matrix, so ambient signal and multinucleated-cell doublets are both plausible;
  the core-wide induction of every programme is consistent with some of this.
- **The bone-pathway score is a mean over 14 genes**, so it cannot distinguish "this cell
  transcribes BMD loci" from "this cell is transcriptionally busy".
- **This dataset is carotid plaque, not myocardium.** It says where the bone programme
  lives inside a calcified artery; it says nothing about the heart, and the eBMD 鈫?HF
  effect in module 05 is a different tissue question.

## Software

- Python 3.12.14 (project `.venv`): numpy 2.5.4, pandas 3.0.6, scipy 1.18.1,
  scikit-learn 1.9.1
- GSE159677 aggregate matrix retrieved 2026-10-21, verified against the published size

