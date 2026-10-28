# 02-mr 鈥?genetic sharing between bone mineral density / osteoporosis and cardiovascular disease

## Purpose

Arbitrate the asymmetry found in `01-nhanes`: at the population level low bone mineral
density behaves more like a shared *prognostic* marker than a shared *prevalent* disease
state. This module asks the genetic question directly 鈥?does a genetically predicted
lower BMD (or genetic liability to osteoporosis) change the risk of coronary artery
disease, myocardial infarction, stroke and heart failure, and does the reverse hold?
Two-sample Mendelian randomisation in both directions is the causal-inference layer (L2)
of the comorbidity study.

## Inputs

Genome-wide summary statistics in `data/raw/gwas/`, obtained from the EBI GWAS Catalog
FTP (`/pub/databases/gwas/summary_statistics/<bin>/<accession>/harmonised/`). Every file
was checked against the byte size and the md5 EBI publishes in `harmonised/md5sum.txt`
(`data/raw/gwas/manifest.csv`, `data/raw/download_manifest.csv`; all six pass).

| Key | Accession | Trait | Role |
|---|---|---|---|
| `eBMD` | GCST90029004 | heel bone mineral density, European (Loh 2018) | exposure |
| `OP` | GCST90018887 | osteoporosis (Sakaue 2021) | exposure |
| `CAD` | GCST90132315 | coronary artery disease (Aragam 2022) | outcome |
| `MI` | GCST90018877 | myocardial infarction (Sakaue 2021) | outcome |
| `STROKE` | GCST90104534 | any stroke (Mishra 2022, GIGASTROKE) | outcome |
| `HF` | GCST90162626 | heart failure (Levin 2022, HERMES) | outcome |

**The files do not share a column schema.** eBMD carries the `hm_*` harmonised columns
plus `n`; the other five carry `rsid / chromosome / base_pair_location / beta /
standard_error / ...` and **no sample-size column**. `scripts/_sumstats.py` normalises
all six onto one schema on read; where `n` is absent the analysis sample size is taken
from the GWAS Catalog study record (`TRAIT_N`) and used **only** for the Steiger test
and instrument R虏, never to derive an effect estimate.

## Methods

1. **Instrument selection** (`scripts/02_extract_instruments.py`). Variants with
   `p < 5e-8`, MAF 1鈥?9 %, valid single-nucleotide alleles, non-missing beta/SE.
2. **LD clumping** (`scripts/02b_ld_clump.py`). A local 1000G phase 3 EUR panel
   (`data/raw/reference/EUR`, PLINK binaries) with `--clump-r2 0.001 --clump-kb 10000`;
   the strongest variant keeps each locus. This replaced an earlier distance-only prune:
   Ensembl's pairwise LD endpoint returns one variant per request, so LD clumping through
   the API was not affordable for 128,923 candidates.
3. **Outcome extraction** (`scripts/03_mr.py`). Each outcome file is streamed once and
   only the instrument SNPs are cached (`data/processed/mr/outcome_snps/`).
4. **Harmonisation.** Outcome betas flipped onto the exposure effect allele; incompatible
   alleles dropped; palindromic A/T and C/G with EAF 0.42鈥?.58 discarded.
5. **Estimators.** IVW (fixed and random), MR-Egger (slope + intercept, SE inflated by the
   residual scale), weighted median and weighted mode (2,000 bootstraps, seeded), simple
   median. Wald ratio by/bx with SE sy/|bx|.
6. **Sensitivity.** Cochran's Q, the MR-Egger intercept, leave-one-out IVW, Steiger.
7. **Pleiotropy correction** (`scripts/04_pleiotropy.py`). Because several pairs carry a
   significant MR-Egger intercept, estimates are re-computed after MR-PRESSO-style
   iterative outlier removal (per-SNP residual z against a Bonferroni threshold), plus a
   leave-one-locus-out sweep over the 10 strongest variants per pair.

## Commands

| Script | Command | Output |
|---|---|---|
| `scripts/00_pipeline_test.py` | `python analysis/02-mr/scripts/00_pipeline_test.py` | end-to-end check on simulated sumstats |
| `scripts/01_download_gwas.py` | `python analysis/02-mr/scripts/01_download_gwas.py --all` | resumable multi-socket download |
| `scripts/01b_verify_manual.py` | `python analysis/02-mr/scripts/01b_verify_manual.py` | verify hand-downloaded files |
| `scripts/02_extract_instruments.py` | `鈥?--datasets eBMD,CAD,MI,STROKE,HF,OP --clump distance --kb 500 --keep-all` | `instruments_*_all.csv` candidates |
| `scripts/02b_ld_clump.py` | `鈥?--datasets OP,STROKE,HF,MI,CAD,eBMD` | `instruments_<KEY>.csv`, `ld_clump_summary.csv` |
| `scripts/03_mr.py` | `python analysis/02-mr/scripts/03_mr.py` | `mr_results_{forward,reverse}.csv` |
| `scripts/04_pleiotropy.py` | `python analysis/02-mr/scripts/04_pleiotropy.py` | `mr_pleiotropy_sensitivity.csv`, `mr_loco_influence.csv` |

## Results as computed

**Instrument sets after LD clumping**

| Dataset | p<5e-8 candidates | LD-independent | min F | 危R虏 |
|---|---|---|---|---|
| eBMD | 128,923 | 600 | 23.4 | 0.182 |
| CAD | 24,015 | 225 | 29.7 | 0.014 |
| MI | 6,722 | 82 | 29.5 | 0.012 |
| STROKE | 1,207 | 44 | 29.7 | 0.002 |
| HF | 917 | 35 | 29.6 | 0.001 |
| OP | 442 | 19 | 29.6 | 0.002 |

危R虏 = 0.182 for eBMD matches published SNP-heritability for heel BMD; every F is far
above the conventional 10 threshold.

**eBMD 鈫?cardiovascular, before and after pleiotropy correction**

| Pair | n | IVW b (95 % CI) | p | Egger intercept p | outliers removed | trimmed b | trimmed p |
|---|---|---|---|---|---|---|---|
| eBMD 鈫?CAD | 591 | 0.033 (0.001, 0.066) | 0.046 | 8e-05 | 19 | 0.021 | **0.13** |
| eBMD 鈫?MI | 592 | 0.026 (鈭?.021, 0.074) | 0.28 | 4.9e-06 | 11 | 0.013 | 0.55 |
| eBMD 鈫?STROKE | 581 | 0.003 (鈭?.024, 0.029) | 0.84 | 0.0053 | 3 | 鈭?.001 | 0.93 |
| eBMD 鈫?HF | 592 | 鈭?.034 (鈭?.059, 鈭?.008) | 0.010 | 0.020 | 4 | **鈭?.040** | **0.0016** |
| eBMD 鈫?OP | 592 | 鈭?.700 (鈭?.760, 鈭?.639) | 7e-114 | 0.25 | 3 | 鈭?.699 | 7e-118 |

**Interpretation.** Only one of the four cardiovascular endpoints survives correction:

- **eBMD 鈫?HF is the single robust finding.** It strengthens after outlier removal
  (鈭?.034 鈫?鈭?.040) and its MR-Egger intercept stops being significant (0.020 鈫?0.055),
  meaning the pleiotropy that inflated the variance is what was removed. Leave-one-locus-
  out over the 10 strongest variants moves the estimate by at most 0.004 (base 鈭?.0399).
  Direction: **higher genetically predicted BMD 鈫?lower heart-failure risk**
  (OR 鈮?0.96 per SD, roughly 4 % per SD).
- **eBMD 鈫?CAD does not survive.** The nominal p = 0.046 disappears after removing 19
  outlier variants (p = 0.13) and the Egger intercept remains significant (p = 0.0024),
  i.e. directional pleiotropy is still present in the trimmed set. This estimate should
  not be reported as a causal effect.
- **eBMD 鈫?MI and eBMD 鈫?STROKE show no effect** at any stage.
- **eBMD 鈫?OP** is the positive control and behaves exactly as expected.

**Reverse direction is not interpretable with these instrument sets.** All CVD 鈫?bone
estimates are fragile: `CAD 鈫?eBMD` flips from p = 0.55 to p = 0.014 only after removing
46 of 205 variants; `MI 鈫?OP` goes from p = 0.037 to p = 0.35 after removing 4 of 78;
`OP 鈫?eBMD` has 17 instruments of which 16 are removed as outliers, leaving the estimate
undefined. With 19鈥?25 instruments whose 危R虏 is 0.001鈥?.014, these sets are too weak for
a confident reverse causal claim, and the outlier procedure itself becomes unstable at
this size.

## Limitations to keep in view

- **Sample overlap.** eBMD is UK Biobank; CAD, MI and HF are largely UK Biobank too.
  Overlap biases two-sample MR toward the observational estimate and makes Steiger
  anti-conservative; Steiger p-values are reported descriptively, not as a filter.
- **Liability-scale exposure.** The `OP` instrument set reflects diagnosed osteoporosis
  (mixing genetic risk with DXA access), and at 19 instruments its estimates are unstable
  under any sensitivity analysis. `OP` results should be read as underpowered, not null.
- **Binary exposure and outcome** are on the log-odds scale; ORs are per log-odds of
  liability, not per SD of BMD.
- **Outlier removal is aggressive at small n.** MR-PRESSO-style trimming removed 16 of 17
  `OP 鈫?eBMD` instruments, which can only mean the method, not the data, is the limit
  there.
- **Long-range LD.** Clumping used a 10 Mb window on a EUR panel; inversions on 8p23 and
  17q21 can still leave correlated instruments. The leave-one-locus-out sweep bounds the
  influence of any single one on the reported pairs (max shift 0.011 for eBMD 鈫?OP).
- **Not done:** LDSC genetic correlation, colocalisation, bidirectional Steiger filtering
  as an inclusion criterion.

## Software

- Python 3.12.14 (project `.venv`): numpy 2.5.4, pandas 3.0.6, scipy 1.18.1,
  matplotlib 3.11.2, statsmodels 0.15.0
- PLINK v1.90b (win64, `data/raw/reference/plink/plink.exe`) with the 1000G phase 3 EUR
  reference (`fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz`, 1.45 GB)
- GWAS Catalog files retrieved 2026-10-21; md5s in `data/raw/gwas/manifest.csv`

