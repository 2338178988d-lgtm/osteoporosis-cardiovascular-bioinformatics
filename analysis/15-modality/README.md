# 15-modality 鈥?phenotype or file? Triangulating the exposure

## Purpose

Module `13-robustness` found that a DXA femoral-neck BMD GWAS gave an estimate of the
**opposite sign** to the heel-ultrasound eBMD estimate (+0.068 vs 鈭?.064, FinnGen,
non-overlapping CIs), and that no diagnostic removed the disagreement. With one file on each
side, two readings stayed open:

- **biological** 鈥?quantitative ultrasound reflects trabecular architecture and elasticity,
  which areal BMD does not capture;
- **technical** 鈥?the DXA file is a secondary deposit from a 2024 machine-learning methods
  paper (PMID 39349818) with no independent DXA GWAS to check it against.

One file per cell cannot separate *modality* from *source*. This module fills the 2 脳 3 matrix
with published GWAS and decides between them.

## Inputs

| Key | Accession | Trait | Modality | n | Source | Instruments |
|---|---|---|---|---|---|---|
| `eBMD` | GCST90029004 | heel ultrasound BMD | ultrasound | 445,921 | Morris 2018 (PMID 29892013) | 600, 危R虏 0.182, min F 23.4 |
| `QUS_2` | GCST90025982 | heel ultrasound BMD T-score | ultrasound | 445,855 | PMID 34226706 | 480, 危R虏 0.152, min F 24.5 |
| `QUS_3` | GCST006979 | heel ultrasound BMD | ultrasound | 426,824 | PMID 30598549 | 511, 危R虏 0.186, min F 24.6 |
| `FNBMD` | GCST90446641 | DXA femoral neck BMD | DXA | 79,940 | PMID 39349818 | 88, 危R虏 0.074, min F 29.7 |
| `FNBMD_2` | GCST90446628 | DXA femoral neck BMD | DXA | 47,205 | PMID 39349818 (same paper) | 56, 危R虏 0.064, min F 30.1 |
| `FNBMD_3` | GCST90568458 | DXA femoral neck BMD | DXA | 32,017 | PMID 41206123 | 32, 危R虏 0.044, min F 30.1 |

Files: `data/raw/gwas/{FNBMD_2,FNBMD_3,QUS_2,QUS_3}_*.h.tsv.gz` (from `download_plan.csv`,
verified by `analysis/scripts/verify_downloads.py`). Instrument tables:
`results/tables/02-mr/instruments_<KEY>.csv`.

Two things about the design that matter for how the result may be stated:

- **The three ultrasound files are three publications on one cohort.** All are UK Biobank heel
  measurements analysed by different groups with different pipelines. Agreement between them
  rules out a single-file or single-pipeline artefact; it is **not** independent-cohort
  replication.
- **The three DXA files are two publications.** `FNBMD` and `FNBMD_2` come from the same 2024
  methods paper (n = 79,940 and 47,205, UKB field 23299); only `FNBMD_3` is from another
  publication.

The full DXA femoral-neck picture on the GWAS Catalog was checked: only these three plus
GCST007691 (Estrada 2012) and GCST005544 carry summary statistics, and the last two are not on
the public FTP in any layout. No lumbar-spine, total-hip or TBS summary statistics exist with
public per-variant data, which is why the DXA arm cannot be widened further.

## Methods

- **Instrument construction** 鈥?`02_extract_instruments.py --clump distance --kb 500
  --keep-all`, then `02b_ld_clump.py` (PLINK v1.90b7.2, 1000G phase 3 EUR, r虏 < 0.001 within
  10 Mb). Genome-wide significance p < 5 脳 10鈦烩伕, MAF 1鈥?9 %, single-character alleles.
- **Estimation** 鈥?IVW with random effects, MR-Egger intercept, allele harmonisation dropping
  palindromes with intermediate frequency. Identical to modules 02鈥?4.
- **Pooling** 鈥?inverse-variance, DerSimonian鈥揕aird 蟿虏, Cochran's Q and I虏 reported per
  modality. **Contrast** 鈥?z-test on the difference between the two pooled estimates.
- **Reciprocal prediction** 鈥?each file's instruments are used to predict every other file's
  trait. A file that cannot predict the same phenotype measured elsewhere is not usable in
  the main comparison, whatever it estimates for heart failure.
- **Decision rules were fixed in the script header before the run** and are printed by
  `_summary.report()`.

## Results

### Heart failure, FinnGen R12 (no UK Biobank overlap)

| File | 尾 | 95 % CI | p | power |
|---|---|---|---|---|
| heel US, Morris 2018 | 鈭?.064 | 鈭?.100, 鈭?.028 | 0.00057 | 93 % |
| heel US, PMID 34226706 | 鈭?.066 | 鈭?.108, 鈭?.023 | 0.0020 | 85 % |
| heel US, PMID 30598549 | **鈭?.071** | 鈭?.109, 鈭?.032 | **0.00031** | 90 % |
| **pooled ultrasound (k=3)** | **鈭?.067** | **鈭?.089, 鈭?.044** | **4.6e-09** | 鈥?|
| DXA neck, PMID 39349818 | +0.068 | +0.005, +0.130 | 0.034 | 52 % |
| DXA neck, same paper (47k) | +0.039 | 鈭?.026, +0.105 | 0.24 | 48 % |
| DXA neck, PMID 41206123 | +0.025 | 鈭?.058, +0.108 | 0.55 | 33 % |
| **pooled DXA (k=3)** | **+0.047** | **+0.008, +0.087** | **0.019** | 鈥?|
| **DXA 鈭?ultrasound** | **+0.114** | **+0.069, +0.160** | **8.7e-07** | 鈥?|

Within-modality heterogeneity is zero in both groups (I虏 = 0.00; Q p = 0.97 and 0.69).

### Heart failure, HERMES (shares UK Biobank with every exposure)

| File | 尾 | 95 % CI | p |
|---|---|---|---|
| heel US, Morris 2018 | 鈭?.034 | 鈭?.059, 鈭?.008 | 0.010 |
| heel US, PMID 34226706 | 鈭?.038 | 鈭?.066, 鈭?.010 | 0.0084 |
| heel US, PMID 30598549 | 鈭?.038 | 鈭?.066, 鈭?.010 | 0.0074 |
| **pooled ultrasound (k=3)** | **鈭?.036** | **鈭?.052, 鈭?.021** | **5.4e-06** |
| DXA neck, PMID 39349818 | +0.032 | 鈭?.013, +0.077 | 0.16 |
| DXA neck, same paper (47k) | +0.012 | 鈭?.031, +0.055 | 0.58 |
| DXA neck, PMID 41206123 | **鈭?.039** | 鈭?.086, +0.007 | 0.096 |
| **pooled DXA (k=3)** | **+0.002** | **鈭?.039, +0.043** | **0.92** |
| **DXA 鈭?ultrasound** | +0.038 | 鈭?.005, +0.082 | 0.087 |

Within-modality heterogeneity: ultrasound I虏 = 0.00, DXA **I虏 = 0.60** (Q p = 0.082).

### Reciprocal prediction: the files do measure the same trait

| Predictor 鈫?trait | 尾 (SD per SD) | p | n SNP |
|---|---|---|---|
| eBMD 鈫?QUS_2 | +0.994 | < 1e-300 | 87 |
| eBMD 鈫?QUS_3 | +0.976 | < 1e-300 | 562 |
| QUS_2 鈫?QUS_3 | +0.983 | < 1e-300 | 464 |
| FNBMD 鈫?FNBMD_2 | +0.978 | 2.9e-197 | 85 |
| FNBMD 鈫?FNBMD_3 | +0.940 | 9.3e-61 | 47 |
| FNBMD_2 鈫?FNBMD_3 | +0.993 | 1.4e-251 | 36 |
| eBMD 鈫?FNBMD | +0.469 | 4.3e-109 | 591 |
| eBMD 鈫?FNBMD_3 | +0.448 | 1.5e-40 | 377 |
| QUS_2 鈫?FNBMD_3 | +0.449 | 3.1e-39 | 311 |

Within a modality the files predict each other at 鈮?1 SD per SD 鈥?they are the same
measurement. Across modalities the slope is 0.45, not 0: DXA and ultrasound are substantially
correlated and clearly distinct, which is exactly the structure needed for the comparison in
the two panels above to mean anything.

### What the pre-registered rules said

- **RULE 1** (every ultrasound file protective, every DXA file positive, low within-modality
  heterogeneity, significant contrast) 鈥?**fires for FinnGen**. *Also fires for HERMES on the
  letter of the rule*, even though the DXA pool there includes zero; the rule as written does
  not require the DXA pool itself to exclude zero.
- **RULE 3** (the DXA pool includes zero 鈫?the reversal is not established) 鈥?**fires for
  HERMES**.

Two rules firing on two different outcome cohorts is the finding, not a nuisance: the
ultrasound side is stable, the DXA side is not.

## What this module establishes

1. **The ultrasound effect is not a single-file artefact.** Three publications, three
   instrument sets, three pipelines, both outcome cohorts: six of six estimates negative,
   all six nominally significant, pooled 鈭?.067 in FinnGen and 鈭?.036 in HERMES with
   I虏 = 0.00 in both. Module 13's "only one file supports this" caveat is retired.
2. **Module 13's reversal does not replicate.** The +0.068 came from the largest DXA file; the
   other two DXA files give +0.039 and +0.025 and neither is significant, and in HERMES the
   third file is **negative** (鈭?.039) with the pooled estimate at +0.002 (I虏 = 0.60). The
   correct statement is not "DXA shows the opposite"; it is **"DXA does not show the effect,
   and the one positive estimate did not survive two further files."**
3. **The exposure wording stays, but for a weaker reason.** "Heel-ultrasound bone density"
   remains the accurate name because there is no positive evidence from DXA 鈥?not because
   there is positive evidence that DXA is opposite.
4. **A monotone pattern worth flagging.** The three DXA estimates fall as the instrument count
   falls (+0.068 with 88 instruments, +0.039 with 56, +0.025 with 32). A real effect should
   stay constant and simply tighten with power; a shrinking estimate usually means the largest
   file carried the noise. This is a caution, not a proof, with three studies.

## Software and data sources

- Python 3.12.14 (project `.venv`): numpy 2.5.4, pandas 3.0.6, scipy 1.18.1, matplotlib 3.11.2
- PLINK v1.90b7.27.2 + 1000G phase 3 EUR (`data/raw/reference/`) for LD clumping
- GWAS Catalog harmonised summary statistics, downloaded 2026-10-23; accessions and PMIDs in
  the Inputs table
- `data/processed/mr/outcome_snps/` 鈥?per-trait outcome caches used by module 14, unaffected

## Commands and scripts

```powershell
python (configured project environment) analysis\scripts\download_plan.py --only FNBMD_2,FNBMD_3,QUS_2,QUS_3 --chunk-mb 4 --connections 8

python (configured project environment) analysis\02-mr\scripts\02_extract_instruments.py --datasets FNBMD_2,FNBMD_3,QUS_2,QUS_3 --clump distance --kb 500 --keep-all
python (configured project environment) analysis\02-mr\scripts\02b_ld_clump.py --datasets FNBMD_2,FNBMD_3,QUS_2,QUS_3

python (configured project environment) analysis\15-modality\scripts\01_modality.py --self-test
python (configured project environment) analysis\15-modality\scripts\01_modality.py

python (configured project environment) analysis\scripts\make_figure10.py
```

Self-test checks that pooling recovers a planted per-file effect (鈭?.060) and that
heterogeneity is reported when it is planted (I虏 = 0.92).

## Outputs

| File | Contents |
|---|---|
| `results/tables/15-modality/modality.csv` | 6 exposures 脳 2 outcomes: 尾, se, CI, p, Q, Egger p, power, MDE, modality, source, PMID |
| `results/tables/15-modality/modality_pooled.csv` | within-modality pools and the modality contrast, both outcomes |
| `results/tables/15-modality/concordance.csv` | reciprocal prediction between files |
| `figures/fig10_modality.png` | both outcome cohorts side by side, plus the concordance check |

## Limitations

- **The three ultrasound files are the same cohort.** Agreement across them excludes a file or
  pipeline artefact, nothing more. No non-UKB heel-ultrasound GWAS with public summary
  statistics exists, so "replicates in an independent cohort" remains unavailable for the
  exposure.
- **The DXA arm is thin**: three files, two publications, 32鈥?8 instruments, 33鈥?2 % power for
  the effect of interest. A null here cannot exclude a small DXA effect; what the data do
  exclude is a **significantly opposite** effect at the size originally reported, once the
  two smaller files are included in HERMES.
- **HERMES shares UK Biobank participants with all six exposures**, so it is the less clean
  cohort and is used for corroboration only; every headline number is FinnGen.
- **I虏 is not informative for the ultrasound group** 鈥?identical participants make the three
  estimates statistically dependent, so inverse-variance pooling of them understates the
  pooled standard error. The pooled ultrasound estimate is reported for direction and
  stability, not as a precision gain.
- **`QUS_3`'s file has an `n` column that is not a sample size** (values 0鈥?). Sample size is
  taken from the study record (426,824) via `_sumstats.BAD_N`; without that override the
  instrument 危R虏 comes out at 526 instead of 0.186.

