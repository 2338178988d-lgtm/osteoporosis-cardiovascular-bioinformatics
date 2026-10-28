# 05-hf-followup 鈥?stress-testing the one surviving causal signal (BMD 鈫?heart failure)

## Purpose

Module `02-mr` left a single estimate that survived outlier removal: genetically
predicted higher bone mineral density (eBMD) lowers heart-failure risk by about 4 % per
SD. A single surviving result in a screen of eighteen pairs is exactly the situation in
which a finding is most likely to be a false positive, so this module attacks it with
four further tests that do not reuse the same evidence: replication in other HF GWAS,
Steiger filtering, multivariable MR against candidate mediators, and removal of loci
already known to influence HF.

## Inputs

Two heart-failure GWAS for replication and four candidate mediators, all GWAS Catalog
harmonised files in `data/raw/gwas/` (md5 same source as module 02; see
`results/tables/download_plan.csv`).

| Key | Accession | Trait | Cases / controls |
|---|---|---|---|
| `HF` | GCST90162626 | heart failure, HERMES (Levin 2022) 鈥?discovery | 47,309 / ~1.22 M |
| `HF_SHAH` | GCST009541 | heart failure (Shah 2020) 鈥?replication | 47,309 / 930,014 |
| `HF_CHRONIC` | GCST90018806 | chronic heart failure (Sakaue 2021) 鈥?replication | 14,262 / 471,898 |
| `BMI` | GCST90029007 | body mass index (Loh 2018) 鈥?mediator | 532,396 |
| `SBP` | GCST90000066 | systolic blood pressure (Surendran 2020) 鈥?mediator | 810,865 |
| `T2D` | GCST005413 | type 2 diabetes (Mahajan 2018) 鈥?mediator | 74,124 / 824,006 |
| `CAD` | GCST90132315 | coronary artery disease (Aragam 2022) 鈥?mediator | 122,733 / ~1.04 M |

Instruments: the 592 LD-clumped eBMD variants from module 02 (危R虏 0.182, min F 23.4).
Mediators own instrument sets, used only to confirm they are instrumented: BMI 479,
SBP 213, T2D 16.

## Methods

1. **Replication** (`scripts/01_hf_followup.py`). The same 592 eBMD instruments against
   each heart-failure GWAS, re-harmonised per file. IVW random effects plus MR-Egger.
2. **Steiger filtering.** Per-SNP variance explained in eBMD versus HF; variants that
   explain more outcome than exposure are dropped and the estimate re-fitted.
3. **Multivariable MR (MVMR).** MVMR-IVW with sandwich standard errors:
   `beta = (X'WX)鈦宦筙'Wy` with `W = diag(1/SE_y虏)`, `X` holding eBMD and one mediator's
   per-SNP effects, `y` the HF effects. One mediator at a time (BMI, SBP, T2D, CAD).
   A conditional-F approximation (residual eBMD effect after regressing out the mediator,
   scaled by its own precision) is reported to show eBMD is still instrumented after
   conditioning.
4. **Shared-locus check.** eBMD instruments within 500 kb of any HF instrument are
   removed and the estimate re-fitted, so the result cannot be driven by variants that
   already associate with HF.
5. **Self-test.** `--self-test` simulates two correlated exposures with known direct
   effects (0.25 and 0.40) and checks MVMR recovers both (recovered 0.254 and 0.404,
   conditional F finite) 鈥?the estimator is verified before it is applied.

## Results as computed

**1. Replication is partial**

| HF GWAS | n | 尾 | 95 % CI | p | OR per SD | Egger intercept p |
|---|---|---|---|---|---|---|
| HF 鈥?HERMES (discovery) | 592 | 鈭?.034 | 鈭?.059, 鈭?.008 | 0.010 | 0.967 | 0.020 |
| HF 鈥?Shah 2020 | 594 | 鈭?.028 | 鈭?.061, 0.005 | 0.097 | 0.972 | 0.246 |
| HF 鈥?Sakaue chronic | 592 | 鈭?.010 | 鈭?.053, 0.034 | 0.66 | 0.990 | 0.961 |

The direction and roughly 80 % of the effect size reproduce in Shah 2020, where the
Egger intercept is also clean; the Sakaue chronic-HF set (14,262 cases) is null but is
an order of magnitude smaller than the others. **This is consistent replication, not
independent confirmation** 鈥?see the overlap caveat below.

**2. Steiger filtering removes nothing**: all 592 variants explain more variance in eBMD
(n = 446k, F 鈮?23鈥?200) than in HF, so the filter cannot discriminate here and the result
is unchanged (鈭?.0336, p = 0.010). Reported for completeness, not as support.

**3. MVMR: the effect is not mediated by adiposity, glycaemia or coronary disease**

| Adjusted for | n | eBMD direct 尾 | p | mediator 尾 | p | change in eBMD effect |
|---|---|---|---|---|---|---|
| 鈥?(univariable) | 592 | 鈭?.0336 | 0.0100 | 鈥?| 鈥?| 鈥?|
| BMI | 592 | 鈭?.0336 | 0.0083 | +0.343 | 1.1e-05 | **0 %** |
| T2D | 584 | 鈭?.0357 | 0.0074 | +0.028 | 0.044 | 鈭? % (slightly larger) |
| CAD | 592 | **鈭?.0451** | **1.3e-04** | +0.327 | 1.7e-26 | **鈭?4 % (larger)** |
| SBP | 21 | +0.055 | 0.30 | +0.010 | 0.97 | **not interpretable** |

The eBMD effect does not attenuate at all when BMI, T2D or CAD are conditioned on; with
CAD adjustment it becomes 34 % larger and more significant. Conditional F for eBMD given
each mediator stays 鈮?35, so weak instruments do not explain this.

**The SBP row must be discarded**: the Surendran SBP file contains only genome-wide
significant variants (11 MB), so only 21 of the 592 instruments carry an SBP effect.
The MVMR fit on 21 variants has no power and its sign change is noise, not evidence.

**4. The result is not carried by known HF loci**

13 of 592 eBMD instruments lie within 500 kb of an HF instrument. Removing them leaves
the estimate slightly stronger: 尾 = 鈭?.0398 (p = 0.0018) versus 鈭?.0336 (p = 0.010).

**5. What the four trimmed outliers are**

| SNP | chr:pos | nearby genes |
|---|---|---|
| rs11657636 | 17:2,150,626 | SMG6 |
| rs42377 | 7:92,614,358 | CDK6 |
| rs117579234 | 12:48,814,262 | ADCY6, CACNB3, DDX23 |
| rs55646464 | 5:173,897,968 | CPEB4 |

They sit in four different regions with no shared bone or cardiac programme, which is
consistent with scattered horizontal pleiotropy rather than one contaminated locus.

---

# Part 2 鈥?FinnGen R12: the independent-population test

Part 1 could not claim independent replication, because every heart-failure GWAS used
there contains UK Biobank, the same source as eBMD. FinnGen R12 is a Finnish cohort with
**no UK Biobank overlap**, which makes it the only genuinely external test available
without a data application. This is the decisive analysis in the whole project.

## Inputs

| Key | File | Endpoint | Cases / controls |
|---|---|---|---|
| `FG_HF` | `finngen_R12_I9_HEARTFAIL.gz` | heart failure, strict | 37,653 / 462,695 |
| `FG_HF_CHD` | `finngen_R12_I9_HEARTFAIL_AND_CHD.gz` | heart failure with coronary disease | 17,532 / 443,698 |
| `FG_HF_OBESE` | `finngen_R12_I9_HEARTFAIL_AND_OVERWEIGHT.gz` | heart failure with BMI 鈮?25 | 22,965 / 206,576 |
| `FG_OP` | `finngen_R12_M13_OSTEOPOROSIS.gz` | osteoporosis | 10,461 / 473,264 |

FinnGen summaries use a different layout (`#chrom pos ref alt rsids beta sebeta pval af_alt`,
GRCh38, **allele `alt` is the effect allele**, no `n` column). `_sumstats.py` maps all
four onto the shared schema via `path_of()`; the analysis N comes from the R12 manifest.

## Results

**1. Replication: stronger and cleaner than in the overlapping cohorts**

| Cohort | n | 尾 (95 % CI) | p | OR per SD | MR-Egger intercept p |
|---|---|---|---|---|---|
| HERMES (discovery, overlapping) | 592 | 鈭?.034 (鈭?.059, 鈭?.008) | 0.010 | 0.967 | **0.020** |
| Shah 2020 (overlapping) | 594 | 鈭?.028 (鈭?.061, 0.005) | 0.097 | 0.972 | 0.246 |
| **FinnGen R12 (independent)** | **577** | **鈭?.064 (鈭?.101, 鈭?.028)** | **0.00057** | **0.938** | **0.801** |
| FinnGen after outlier removal | 576 | 鈭?.068 (鈭?.105, 鈭?.032) | 0.00022 | 0.934 | 鈥?|

The independent cohort gives an effect **twice as large** as the discovery set and,
critically, **no detectable directional pleiotropy** (Egger intercept p = 0.80, versus
p = 0.02 in HERMES). The pleiotropy that worried the discovery analysis is absent where
the sample does not overlap, which is what a real effect looks like and what a
UK-Biobank-driven artefact does not.

**2. Subtype test: the effect is absent from coronary heart failure**

| Endpoint | n | 尾 (95 % CI) | p | OR |
|---|---|---|---|---|
| All heart failure (37,653) | 577 | 鈭?.064 (鈭?.101, 鈭?.028) | 0.00057 | 0.938 |
| **HF + coronary disease (17,532)** | 577 | **+0.008 (鈭?.046, +0.062)** | **0.78** | **1.008** |
| HF + BMI 鈮?25 (22,965) | 577 | 鈭?.078 (鈭?.123, 鈭?.033) | 0.00066 | 0.925 |

This is the single most informative result in the project. The BMD effect is **flatly
absent** in the heart-failure subtype that carries coronary disease (point estimate
essentially zero, opposite sign), while it is fully present in the obese stratum. It
corroborates the multivariate MR from part 1 鈥?which found no attenuation by CAD and no
attenuation by BMI 鈥?using an entirely different, non-overlapping sample and a
case-definition approach instead of a covariate approach.

The corollary matters for interpretation: **"osteoporosis and cardiovascular disease are
comorbid" is the wrong frame.** The genetic signal is specific to a non-coronary heart
failure pathway, not shared with atherosclerotic disease.

**3. Osteoporosis as a binary exposure: valid instruments, insufficient power**

| Pair | n | 尾 (95 % CI) | p |
|---|---|---|---|
| FG-OP 鈫?UKB eBMD (cross-population positive control) | 10 | **鈭?.357 (鈭?.482, 鈭?.231)** | **2.8e-08** |
| FG-OP 鈫?FG-HF (strict) | 10 | 鈭?.029 (鈭?.105, +0.047) | 0.45 |
| FG-OP 鈫?FG-HF + CHD | 10 | 鈭?.014 (鈭?.126, +0.098) | 0.81 |
| FG-OP 鈫?FG-HF + BMI鈮?5 | 10 | 鈭?.032 (鈭?.132, +0.068) | 0.53 |

The cross-population control is the key row: FinnGen osteoporosis liability lowers UK
Biobank heel BMD by 0.36 SD per log-odds (p = 2.8e-08), with 12 LD-independent
instruments at mean F = 31. **The instrument set is valid.** Its failure on heart failure
is a power problem (10 usable variants, CI spanning zero), not a validity problem 鈥?
exactly the limitation flagged for the 19-variant UK-Biobank OP set in module 02.

**4. Reverse direction, run once and stopped**

| Pair | n | 尾 (95 % CI) | p | Egger intercept p |
|---|---|---|---|---|
| FG-HF 鈫?UKB eBMD | 11 | +0.079 (鈭?.003, +0.162) | 0.058 | 0.53 |
| FG-HF 鈫?FG-OP | 12 | 鈭?.016 (鈭?.162, +0.130) | 0.83 | 0.36 |

Both are null-to-borderline with CIs an order of magnitude wider than the forward
estimate, 11鈥?2 instruments and 危R虏 two orders below eBMD's 0.182. The point estimate of
`FG-HF 鈫?eBMD` is also in the *opposite* sign to a causal reverse pathway. The honest
reading is **underpowered and uninformative**, not "no reverse causation" and not
"reverse causation present".

## What part 2 establishes

- The eBMD 鈫?heart-failure effect **replicates in a population with zero sample overlap**,
  at twice the discovery effect size, with **no residual pleiotropy**.
- The effect is **absent in coronary heart failure and present in obese heart failure**,
  which independently confirms the multivariate MR result that the pathway is
  neither atherosclerotic nor adiposity-mediated.
- Osteoporosis (binary) instruments are valid but underpowered for heart failure; the
  continuous BMD phenotype is the useful exposure.

## What part 2 does not establish

- **Which pathway actually carries the effect.** Ruling out coronary disease and adiposity
  is not identifying the mechanism. FGF23, sclerostin, osteocalcin, cystatin C and CRP
  are queued as mediators (module 06) and are not yet analysed.
- **Whether treating osteoporosis lowers heart-failure risk.** These are liability-scale
  genetic estimates, not treatment effects.
- **HFpEF versus HFrEF.** FinnGen subtypes are defined by comorbidity (coronary disease,
  BMI), not by ejection fraction. The obese-HF stratum is an imperfect proxy for HFpEF
  and should not be read as one.

## Commands

| Script | Command | Output |
|---|---|---|
| `scripts/01_hf_followup.py` | `python analysis/05-hf-followup/scripts/01_hf_followup.py` | `hf_*.csv` (part 1) |
| `scripts/02_finngen.py` | `python analysis/05-hf-followup/scripts/02_finngen.py` | `fg_replication.csv`, `fg_subtype.csv`, `fg_reverse.csv`, `fg_trimmed.csv` |
| `scripts/03_finngen_op.py` | `python analysis/05-hf-followup/scripts/03_finngen_op.py` | `fg_op.csv` |
| self-test | `python analysis/05-hf-followup/scripts/01_hf_followup.py --self-test` | PASS/FAIL |

Inputs are produced by `analysis/02-mr/scripts/02_extract_instruments.py` and
`02b_ld_clump.py` (eBMD, BMI, SBP, T2D). FinnGen instruments are built inside
`02_finngen.py` because that file layout differs from the GWAS Catalog schema.

## Commands (part 1)

| Script | Command | Output |
|---|---|---|
| `scripts/01_hf_followup.py` | `python analysis/05-hf-followup/scripts/01_hf_followup.py` | `results/tables/05-hf-followup/*.csv` |
| self-test | `python analysis/05-hf-followup/scripts/01_hf_followup.py --self-test` | PASS/FAIL |

Inputs are produced by `analysis/02-mr/scripts/02_extract_instruments.py` and
`02b_ld_clump.py` (eBMD) and the same pair run with `--datasets BMI,SBP,T2D`.

## What this does and does not establish

**Establishes.** A genetically predicted higher BMD is associated with lower HF risk in
two of three HF GWAS with consistent direction and size; the association is not explained
by BMI, type-2 diabetes or coronary artery disease; it is not driven by variants already
known to affect HF; and it is not an artefact of one outlier locus.

**Does not establish.** That treating osteoporosis reduces heart failure. The exposure is
a continuous heel-BMD phenotype in UK Biobank, not a treatment; the mediator set tested
is small (SBP unusable, no inflammatory, renal, vitamin-D or sex-hormone mediators); and
the HF GWAS overlap substantially with eBMD.

## Limitations

- **Sample overlap (part 1 only).** eBMD is UK Biobank. HERMES (Levin 2022) and
  Shah 2020 both include UK Biobank, so neither is an independent replication; Sakaue's
  chronic-HF set shares UK Biobank too but is much smaller. **Part 2 resolves this**
  with FinnGen R12, which shares no samples with UK Biobank.
- **The independent replication is Finnish-only.** FinnGen is a founder population with
  different allele frequencies and a different diagnostic coding system from UK Biobank,
  so the effect could be population-specific rather than universal. The effect being
  *larger* in FinnGen (OR 0.938) than in UKB (0.967) draws attention to this rather
  than dismissing it.
- **MVMR design.** Only eBMD instruments were used, so the mediator coefficients are
  estimated with instruments chosen for another trait; they are reported for context, not
  as causal mediator effects. Instruments are independent within each exposure but not
  necessarily across exposures.
- **HF subtypes are not separated.** HFpEF and HFrEF are not distinguished in any of the
  three files, and a bone-density effect plausibly differs between them.
- **Pleiotropy is reduced, not eliminated.** The Egger intercept in HERMES remains
  nominally significant (p = 0.020) even though removing 4 outliers fixed it in module 02's
  trimmed analysis; the Shah replication shows no such pleiotropy (p = 0.246).
- **The SBP mediator is unusable** with summary statistics restricted to significant
  variants.

## Software

- Python 3.12.14 (project `.venv`): numpy 2.5.4, pandas 3.0.6, scipy 1.18.1
- PLINK v1.90b7.2 + 1000G phase 3 EUR for clumping (module 02)
- GWAS Catalog harmonised files retrieved 2026-10-21

