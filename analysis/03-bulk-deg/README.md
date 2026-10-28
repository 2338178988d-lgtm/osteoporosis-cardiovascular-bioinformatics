# 03-bulk-deg 鈥?two-tissue differential expression (bone side 脳 vascular side)

## Purpose

Layer L3 of the comorbidity study. The design constraint established in
`results/reports/鍏叡鏁版嵁璧勬簮娓呭崟涓庣爺绌舵灦鏋刜v1.md` is that **no public dataset contains both
bone and vascular transcriptomes from the same people**. This module therefore takes the
design the constraint allows: find differentially expressed genes separately on each
side, then ask whether the two sides converge 鈥?first at the gene level, then at the
pathway level.

**Outcome: they do not converge.** The negative result is the module's actual finding and
is quantified below rather than hidden behind a threshold.

## Inputs

13 GEO series matrices in `data/raw/geo/<GSE>/` (93 MB total, verified with
`analysis/scripts/verify_downloads.py`), plus two supplementary files
(GSE159677 scRNA-seq tar, GSE221911 TPM xlsx) used by other modules.

| GSE | Species | Design | n | Role |
|---|---|---|---|---|
| GSE56815 | human | circulating monocytes, low vs high hip BMD | 80 | bone discovery |
| GSE56814 | human | same cohort, HuEx-1_0-st | 73 | platform replicate |
| GSE7158 | human | monocytes, extremely low vs high peak bone mass | 26 | bone validation |
| GSE7429 | human | circulating B cells, low vs high BMD | 20 | bone immune subset |
| GSE13850 | human | circulating B cells, low vs high BMD (smokers) | 40 | bone immune subset |
| GSE35956 | human | femoral-head hMSC, primary osteoporosis vs matched | 10 | bone tissue |
| GSE100927 | human | carotid/femoral/infra-popliteal, atherosclerotic vs control | 104 | vascular discovery |
| GSE43292 | human | carotid atheroma vs intact intima, paired per patient | 64 | vascular replicate |
| GSE28829 | human | advanced vs early carotid plaque | 29 | vascular progression |
| GSE37558 | human | calcifying vascular cells vs osteoblasts (in vitro) | 32 | mechanism |
| GSE11917 | human | coronary SMC calcification time course | 105 | mechanism (not analysed: gradient design) |
| GSE80614 | human | hMSC osteogenic vs adipogenic time course | 69 | mechanism (not analysed: gradient design) |
| GSE51495 | **baboon** | PBMC vs cortical bone, same animals | 30 | **excluded: not human** |

> **GSE51495 correction.** The v1 data report described GSE51495 as "the same people's
> PBMC + cortical bone". Reading the downloaded series matrix shows all 30 samples are
> `Baboon cortical bone` / `Baboon Peripheral blood mononuclear cells`. It is therefore
> **not usable as evidence that human blood substitutes for bone**, and the study design
> loses its only human blood鈥揵one bridge. The correction is recorded in
> `config/datasets_meta.csv` and in the data report.

## Methods

1. **Series-matrix parsing** (`scripts/01_series_matrix.R`). Metadata lines plus the
   `!series_matrix_table` are read directly. The assay is log2-transformed **only if the
   99th percentile is 鈮?30** (so already-log2 data are not transformed twice), then
   quantile-normalised (`limma::normalizeBetweenArrays`). Quote characters are stripped
   from probe IDs 鈥?the matrix header is quoted and `quote=""` otherwise keeps them,
   which silently breaks every probe鈫抔ene lookup.
2. **Probe 鈫?gene** (`scripts/00_build_symbols.py`). Platform annotation from
   `annot/GPLxxx.annot.gz` when published (GPL96, GPL570, GPL6104, GPL6244, GPL6947,
   GPL17077), otherwise from the platform table via `acc.cgi?...&view=data`
   (GPL5175). Duplicate symbols keep the highest-mean probe. All 13 datasets map
   (13,299鈥?0,558 genes each).
3. **Differential expression** (`scripts/02_deg_limma.R`). Groups come from
   `config/datasets.csv` (per-dataset regexes on sample title + characteristics).
   `lmFit` + `eBayes(trend = TRUE)`, `topTable(coef = "groupcase")`. GSE43292 is a paired
   design (plaque and intact intima from the same patient) and is fitted with
   `duplicateCorrelation` + `block = patient`.
4. **Intersection** (`scripts/03_intersect.py`). Direction-consistent overlap between the
   bone side and the vascular side at six thresholds, each tested against the
   hypergeometric expectation. Within-side replication is computed the same way.
5. **Pathway convergence** (`scripts/04_pathway_convergence.R`). ORA against 50 MSigDB
   Hallmark sets, restricted to the genes measurable on both sides (10,014 genes), then
   the number of shared enriched pathways is compared with a permutation null over set
   names.
6. **Self-tests.** Every script has a `--selftest` (or `--self-test`) mode on synthetic
   data: quantile normalisation equalises sample medians, log2 decision fires on
   linear input, 37/40 planted DEGs land in the top 10 %, and the permutation test gives
   p = 0.001 on a genuine pathway overlap and p = 1.0 on a disjoint one. All pass.

## Results as computed

**Differential expression (FDR < 0.05, |log2FC| > 1)**

| Dataset | Side | genes tested | FDR<0.05 | up | down | max |log2FC| |
|---|---|---|---|---|---|---|
| GSE100927 | vascular | 30,558 | 14,056 | 388 | 152 | 3.85 |
| GSE43292 | vascular | 20,358 | 8,125 | 61 | 57 | 2.43 |
| GSE28829 | vascular | 22,189 | 2,229 | 287 | 41 | 4.14 |
| GSE37558 | mechanism | 19,613 | 4,872 | 200 | 144 | 4.70 |
| GSE35956 | bone | 22,189 | 401 | 163 | 49 | 3.97 |
| GSE56815 | bone | 13,299 | 432 | 0 | 0 | **0.68** |
| GSE56814 | bone | 17,324 | 0 | 0 | 0 | **0.74** |
| GSE7158 | bone | 22,189 | 0 | 0 | 0 | **1.77** |
| GSE7429 | bone | 13,299 | 0 | 0 | 0 | 0.77 |
| GSE13850 | bone | 13,299 | 1,139 | 2 | 0 | 2.44 |

**Within-side replication (FDR<0.05, |log2FC|>0.5; observed vs chance)**

| Pair | observed | expected | fold | p |
|---|---|---|---|---|
| GSE100927 vs GSE43292 | 767 | 204.2 | 3.76 | 1.7e-312 |
| GSE100927 vs GSE28829 | 627 | 175.1 | 3.58 | 4.5e-233 |
| GSE43292 vs GSE28829 | 557 | 85.0 | 6.55 | ~0 |
| GSE56815 vs GSE56814 | 0 | 鈥?| 鈥?| (2 vs 0 DEGs) |
| GSE56815 vs GSE7158 | 0 | 鈥?| 鈥?| (2 vs 0 DEGs) |
| GSE13850 vs GSE35956 | 13 | 2.84 | 4.58 | 4.7e-06 |

**Vascular replication is strong and consistent. Bone-side blood studies are not: the
same-cohort platform replicate (GSE56815 vs GSE56814) shares zero DEGs, and five of six
bone datasets produce fewer than 3 genes past |log2FC| > 1.**

**Bone 脳 vascular intersection (direction-consistent, hypergeometric null)**

| Threshold | bone | vascular | observed | expected | fold | p |
|---|---|---|---|---|---|---|
| FDR<0.05, |FC|>1 | 0 | 502 | 0 | 鈥?| 鈥?| 鈥?|
| FDR<0.05, |FC|>0.5 | 2 | 2,229 | 0 | 0.18 | 0 | 1.0 |
| FDR<0.05 | 432 | 10,510 | 112 | 188.2 | 0.60 | 1.0 |
| FDR<0.10 | 906 | 11,374 | 251 | 427.1 | 0.59 | 1.0 |
| P<0.05 | 2,950 | 11,349 | 798 | 1,387.6 | 0.58 | 1.0 |
| P<0.01 | 1,135 | 9,743 | 280 | 458.3 | 0.61 | 1.0 |

**At every threshold the overlap is below chance, never above.** The 56 genes that do
pass FDR<0.05 with |FC|>0.5 in both directions (`intersection_genes.csv`) are therefore
not evidence of sharing.

**Pathway convergence (Hallmark, FDR<0.05)**

| Side | enriched pathways |
|---|---|
| vascular | 24 |
| mechanism (GSE37558) | 7 |
| bone | 4 |

Shared bone 鈭?vascular = **2** (HALLMARK_MYC_TARGETS_V1, HALLMARK_E2F_TARGETS) against
**1.92 expected**; permutation **p = 0.66**. Bone 鈭?mechanism = 1 (p = 0.45);
vascular 鈭?mechanism = 5 (p = 0.18). The two shared pathways are generic proliferation
programmes, not a bone鈥搗ascular axis.

## Commands

| Script | Command | Output |
|---|---|---|
| `scripts/00_build_symbols.py` | `python analysis/03-bulk-deg/scripts/00_build_symbols.py` | `data/processed/geo/_platform/<GSE>_symbols.tsv` |
| `scripts/01_series_matrix.R` | `Rscript analysis/03-bulk-deg/scripts/01_series_matrix.R` | `data/processed/geo/<GSE>/{expression,phenotype,qc}.csv` |
| `scripts/02_deg_limma.R` | `Rscript analysis/03-bulk-deg/scripts/02_deg_limma.R` | `results/tables/03-bulk-deg/deg_<GSE>.csv`, `deg_summary.csv` |
| `scripts/03_intersect.py` | `python analysis/03-bulk-deg/scripts/03_intersect.py` | `intersection_grid.csv`, `intersection_genes.csv`, `within_side_replication.csv` |
| `scripts/04_pathway_convergence.R` | `Rscript analysis/03-bulk-deg/scripts/04_pathway_convergence.R` | `pathway_ora_*`, `pathway_summary.csv`, `pathway_convergence.csv` |
| self-tests | `--selftest` on each | PASS/FAIL |

## Limitations

- **The bone side is blood, not bone.** Five of six bone datasets are circulating
  monocytes or B cells; only GSE35956 (n=10) is bone-derived, and it is underpowered.
  The surrogate assumption has **no human blood鈥揵one pairing to support it** (GSE51495 is
  baboon), so the bone side of every comparison is an assumption, not a measurement.
- **Weak bone-side signal is not purely a sample-size problem.** GSE56814 (n=73) and
  GSE7158 (n=26) produce zero genes past |log2FC|>1 while GSE100927 (n=104) produces 540;
  the difference tracks tissue and phenotype contrast (atherosclerotic plaque vs intact
  artery is a large biological difference; low vs high BMD monocytes is not).
- **GSE7158's `source_name` is constant** ("circulating monocytes in low-pbm phase"), which
  is why grouping uses title + characteristics only.
- **GSE43292 pairing** is modelled with `duplicateCorrelation`; a mixed-effects model
  would be equivalent but was not cross-checked.
- **GSE11917 and GSE80614 are gradient/time-course designs** and are still not analysed;
  they need a polynomial or spline model rather than a two-group contrast.
- **Cross-platform comparison.** The sides share only 10,014 genes across five different
  platforms; gene-level intersection is therefore on a restricted, platform-biased
  universe.
- **No WGCNA, no single-cell localisation, no GTEx baseline** was run, although the data
  (GSE159677) are downloaded.

## Software

- R 4.6.1 (`Rscript (configured environment)`): limma 3.68.5, clusterProfiler 4.20.0,
  msigdbr 26.1.1, org.Hs.eg.db 3.23.1, ReactomePA 1.56.0, GEOquery 2.80.0, WGCNA 1.74,
  survey 4.5, survival 3.8.6
- Python 3.12.14 (project `.venv`): numpy 2.5.4, pandas 3.0.6, scipy 1.18.1
- GEO series matrices retrieved 2026-10-21; md5s in `data/raw/download_manifest.csv`


