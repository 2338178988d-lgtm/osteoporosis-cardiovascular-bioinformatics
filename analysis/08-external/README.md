# 08-external 鈥?a third population, and ejection fraction

## Purpose

Module `05-hf-followup` established the eBMD 鈫?heart-failure effect and replicated it in
FinnGen. Two questions were left open, and both are addressed here.

**1. Is the effect a Finnish peculiarity?** The FinnGen effect was *larger* than in UK
Biobank (OR 0.938 vs 0.967), which invites the objection that a founder population with a
different diagnostic coding system produced a different answer. A third population with no
UK Biobank and no Finnish overlap settles it 鈥?or fails to.

**2. Does ejection fraction matter?** Module 05's subtype test was by *comorbidity*
(coronary disease, obesity) because FinnGen publishes no ejection-fraction split. That
test found the effect absent from coronary heart failure. An ejection-fraction split asks a
different and more mechanistic question.

## Inputs

| Key | Accession | Trait | Cases / controls | Population |
|---|---|---|---|---|
| `FG_EA_HF` | GCST90018586 | chronic heart failure (Sakaue 2021) | 178,726 total | East Asian (BioBank Japan) |
| `FG_EA_OP` | GCST90018667 | osteoporosis (Sakaue 2021) | 178,726 total | East Asian |
| `HFPEUR` | GCST90726619 | non-ischaemic HFpEF (Henry 2025) | 3,590 / 752,599 | European |
| `HFREUR` | GCST90726620 | non-ischaemic HFrEF (Henry 2025) | 4,975 / 935,830 | European |
| `HFPEA` | GCST90668011 | HFpEF (Enzan 2025) | 鈥?| East Asian |
| `HFREA` | GCST90668010 | HFrEF (Enzan 2025) | 鈥?| East Asian |

Exposure: the 600 LD-clumped eBMD instruments. Comparators: `HF` (HERMES) and `FG_HF`
(FinnGen) from module 05. All six files are now present and gzip-verified; sizes and md5
in `data/raw/download_manifest.csv`.

## Results

### Third population: BioBank Japan is null, and underpowered

| Cohort | n | 尾 (95 % CI) | p | OR per SD | Egger intercept p | power to detect 鈭?.064 |
|---|---|---|---|---|---|---|
| UKB / HERMES (overlapping) | 592 | 鈭?.034 (鈭?.059, 鈭?.008) | 0.010 | 0.967 | 0.020 | 鈥?|
| FinnGen R12 (Finnish) | 577 | 鈭?.064 (鈭?.100, 鈭?.028) | 0.00057 | 0.938 | 0.801 | ~99 % |
| **BioBank Japan (East Asian)** | **468** | **鈭?.027 (鈭?.095, +0.041)** | **0.44** | **0.973** | **0.0035** | **46 %** |

Cochran's Q across the three cohorts is 2.01 (df 2, p = 0.37): no evidence of
heterogeneity, so the null is not a contradiction of the Finnish result. But it is not a
replication either. The BioBank Japan estimate sits on the same side of zero as the other
two and its confidence interval includes the FinnGen effect; with 468 instruments and a
standard error of 0.035 the study had 46 % power to detect 鈭?.064. Calling this a failed
replication would be as wrong as calling it a confirmation.

Two further reasons not to over-read it: the FinnGen instrument definitions come from a
European eBMD GWAS applied to an East Asian outcome (allele frequencies and LD differ),
and the BioBank Japan file shows directional pleiotropy on its own (MR-Egger intercept
p = 0.0035), unlike FinnGen (p = 0.80).

**The reverse East Asian test cannot be run.** Genome-wide significant BioBank Japan heart
failure variants number 10, of which 1 survives LD clumping (min F = 30.9) 鈥?below the
10-instrument floor the module requires, so `ea_reverse.csv` is not produced. The reverse
direction stays unanswered, as it was in module 05.

### Ejection fraction: European HFpEF is now included, and it is still a null with a suggestive shape

| Subset | n | 尾 (95 % CI) | p | SE | power to detect the all-HF effect |
|---|---|---|---|---|---|
| **HFpEF, European** (3,590 cases) | 573 | **鈭?.056 (鈭?.158, +0.046)** | 0.28 | 0.052 | **23 %** |
| HFrEF, European (4,975 cases) | 574 | +0.005 (鈭?.086, +0.096) | 0.92 | 0.047 | **28 %** |
| HFpEF, East Asian | 462 | +0.004 (鈭?.081, +0.089) | 0.93 | 0.043 | **31 %** |
| HFrEF, East Asian | 462 | 鈭?.022 (鈭?.126, +0.082) | 0.67 | 0.053 | **23 %** |
| 鈥?all HF, FinnGen (reference) | 577 | 鈭?.064 (鈭?.101, 鈭?.028) | 0.0006 | 0.019 | ~99 % |

HFpEF versus HFrEF: European difference 鈭?.061 (p 鈮?0.38), East Asian +0.026 (p 鈮?0.70).

**The European pair has the shape the rest of the project predicts** 鈥?HFpEF carries the
full all-HF effect (鈭?.056 against 鈭?.064) while HFrEF is flat (+0.005) 鈥?and that is
consistent with module 05's finding that the effect disappears from heart failure defined
with coronary disease and survives in obesity-defined disease. It is nonetheless **not
evidence**: neither estimate is individually significant, the contrast does not approach
significance, and the observed power is 23鈥?8 %. The same shape from a 20,000-case HFpEF
GWAS would be a mechanism; from 3,590 cases it is a hypothesis.

A second limitation stands: Henry 2025's HFpEF and HFrEF are **non-ischaemic** by
construction, while the East Asian files are not, so the two populations are not testing
the same phenotype.

## What this module establishes

It closes the two questions it was built for, and neither closes in the paper's favour:

- The ejection-fraction test **cannot be answered** at current GWAS sizes: 3,590鈥?,975
  cases gives 23鈥?1 % power for the effect size in question, and the one suggestive
  pattern (European HFpEF) is indistinguishable from noise.
- The third-population test **has been run and is inconclusive**: BioBank Japan shows no
  effect (p = 0.44) with 46 % power, no heterogeneity across the three cohorts (p = 0.37),
  and no usable instruments for the reverse direction.

## What would make these answerable

| Question | Requirement | Status |
|---|---|---|
| Ejection fraction | an HFpEF/HFrEF GWAS with 鈮?20,000 cases each | not available publicly at that size |
| Third population | a non-European, non-Finnish HF GWAS with 鈮?2脳 the cases of GCST90018586, plus a matched-ancestry exposure | GCST90018586 is the largest chronic-HF file in BioBank Japan |
| Both | pooled multi-ancestry HF GWAS with subtype annotation | GCST90726617 exists (139,533 European cases) but is not stratified by EF |

## Commands

| Script | Command | Output |
|---|---|---|
| `scripts/01_external_validation.py` | `python analysis/08-external/scripts/01_external_validation.py` | `ea_replication.csv`, `ea_heterogeneity.csv`, `ef_subtype.csv`, `ef_contrast.csv` |
| self-test | `... --self-test` | PASS/FAIL |

Self-test checks the cohort-heterogeneity statistic (identical cohorts Q = 0, p = 1;
divergent cohorts Q = 165, p = 1e-36).

The script skips any trait whose file is absent rather than failing, so it can be re-run as
downloads complete. It was written for exactly that reason: the first run (2026-10-22,
before the two large files finished) produced only the overlapping cohorts and three
ejection-fraction subsets.

## Limitations

- **Power, not absence.** See above; this is the dominant limitation and it invalidates any
  claim about ejection fraction or about the East Asian cohort.
- **Population mismatch in the subtype files.** European HFpEF/HFrEF are non-ischaemic by
  construction; the East Asian files are not, so the two populations are not testing the
  same phenotype.
- **No East Asian instrument set was built.** Ideally the exposure would be East Asian
  eBMD for an East Asian outcome; no large East Asian heel-BMD GWAS with summary
  statistics was found, so a European instrument set is applied to an East Asian outcome.
  That is defensible for a locus with consistent cross-ancestry effects but is an
  approximation, and it is the leading alternative explanation for the BioBank Japan null.
- **The reverse East Asian test remains unrun**: 10 significant BioBank Japan heart-failure
  variants collapse to 1 independent instrument, below the module's floor of 10.
- **BioBank Japan's own pleiotropy.** Its MR-Egger intercept is significant (p = 0.0035)
  while FinnGen's is not (p = 0.80), so part of the difference between the two may be
  instrument behaviour rather than population.

## Software

- Python 3.12.14 (project `.venv`): numpy 2.5.4, pandas 3.0.6, scipy 1.18.1
- Files retrieved 2026-10-22 from the EBI GWAS Catalog; sizes and md5 in
  `results/tables/download_plan.csv`

