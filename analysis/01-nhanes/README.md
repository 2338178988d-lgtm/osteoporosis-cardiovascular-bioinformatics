# 01-nhanes 鈥?bone mineral density 脳 cardiovascular disease in US adults

## Purpose

Test, in a nationally representative sample, whether lower femoral neck bone
mineral density (BMD) is associated with (i) prevalent cardiovascular disease
(CVD) and (ii) subsequent all-cause and cardiovascular mortality. This module is
the epidemiology layer (L1) of the osteoporosis鈥揅VD comorbidity study; it
establishes whether a comorbidity signal exists at all before any transcriptomic
work is commissioned.

## Inputs

Raw NHANES files downloaded by `scripts/01_download.py` into
`data/raw/nhanes/<cycle>/` (immutable; URL, bytes, SHA-256 and row count for each
file are recorded in `data/raw/nhanes/manifest.csv`, 125 files):

| File group | Role |
|---|---|
| `DEMO_*` | age, sex, race/ethnicity, education, poverty-income ratio, MEC weight, PSU, strata |
| `DXXFEM_*` | femoral neck / total femur / trochanter BMD 鈥?primary exposure (cycles D, E, F, H, J) |
| `DXXSPN_*` | lumbar spine BMD 鈥?sensitivity exposure |
| `DXX_*` | whole-body DXA file that replaces the split femur/spine files in cycles G and I; contains no proximal-femur BMD, so G and I cannot contribute a T-score |
| `OSQ_*` | self-reported osteoporosis (absent in 2011鈥?016) |
| `MCQ_*` | self-reported CHF, CHD, angina, MI, stroke 鈥?CVD outcome |
| `BPQ_*`, `BPX_*` | hypertension diagnosis/medication and measured blood pressure |
| `DIQ_*`, `GHB_*` | diabetes diagnosis and HbA1c |
| `SMQ_*`, `BMX_*`, `TCHOL_*`, `HDL_*`, `BIOPRO_*`, `RHQ_*` | covariates |
| `mortality/NHANES_<cycle>_MORT_2019_PUBLIC.dat` | linked mortality (fixed width), follow-up to 2019-12-31 |

Mortality file layout was not guessed: it was read from the CDC public-use 2019
LMF data dictionary (`literature/nhanes-lmf/public-use-linked-mortality-files-data-dictionary.pdf`)
and then confirmed empirically against the downloaded files (column-value
distributions are consistent with MORTSTAT 0/1, UCOD_LEADING 001鈥?10,
PERMTH 0鈥?74).

## Methods

1. **Harmonisation** (`scripts/02_harmonize.py`). Five cycles with proximal-femur
   DXA (2005鈥?6, 2007鈥?8, 2009鈥?0, 2013鈥?4, 2017鈥?8) are stacked and merged by
   `SEQN`. 2011鈥?2 and 2015鈥?6 are excluded because their DXA file carries no
   femoral neck measurement.
2. **T-score reference.** Sex-specific young-adult reference (age 20鈥?9) was
   computed *within these cycles* rather than taken from published constants:
   femoral neck BMD 0.973 卤 0.144 g/cm虏 (men, n=1,149) and 0.890 卤 0.132 g/cm虏
   (women, n=1,062). T-score = (BMD 鈭?mean)/SD. WHO categories: 鈮?鈭?.5
   osteoporosis, 鈭?.5 to 鈭?.0 osteopenia, > 鈭?.0 normal.
3. **Outcomes.** Prevalent CVD = any self-reported physician diagnosis of CHF,
   CHD, angina, MI or stroke. Mortality outcomes from the linked file:
   all-cause death, and CVD death (UCOD_LEADING 001 heart disease or 005
   cerebrovascular disease).
4. **Analytic sample.** Age 鈮?40, non-missing femoral neck BMD, non-missing CVD
   status, mortality-eligible (`ELIGSTAT` = 1): n = 14,260 (flow in
   `table0_sample_flow.csv`).
5. **Survey design** (`scripts/03_survey_analysis.R`). `svydesign(ids = ~psu,
   strata = ~cycle:strata, weights = WTMEC2YR/5, nest = TRUE)`; all estimates and
   models are design-based; sex and age subgroups use domain analysis on the full
   design, never pre-filtered data frames.
6. **Models.** Logistic (quasi-binomial) regression for prevalent CVD with four
   nested adjustment sets: M1 crude; M2 + age, sex, race, cycle; M3 + education,
   PIR, BMI, smoking, diabetes, hypertension, total cholesterol, HDL, eGFR;
   M4 + albumin, calcium, phosphate, alkaline phosphatase. Cox models
   (`svycoxph`) for mortality with the same escalation plus prevalent CVD.
   Continuous exposures: 鈭扵-score (per 1 SD lower) and BMD per 0.1 g/cm虏 lower.
   Linear trend tested across WHO categories.
7. **Dose-response** (`scripts/04_rcs.R`): natural cubic spline (3 df) of
   T-score, OR referenced to T = 0, contrasts computed from the design-based
   covariance matrix.
8. **Survival curves** (`scripts/05_km.R`): weighted Kaplan-Meier from
   Horvitz鈥揟hompson risk sets with weights normalised to mean 1 (raw survey
   weights collapse the Greenwood variance and produced zero-width intervals).
   Descriptive only; inference comes from the Cox models.

## Software and data sources

- Python 3.12.14 (project `.venv`): numpy 2.5.4, pandas 3.0.6, scipy 1.18.1,
  matplotlib 3.11.2, statsmodels 0.15.0, scikit-learn 1.9.1, Pillow 13.0.0
- R 4.6.1 (`Rscript (configured environment)`): survey 4.5, survival 3.8-6,
  splines (base), data.table 1.18.6.1
- NHANES: CDC/NCHS continuous NHANES public release, files retrieved 2026-10-21
  from `https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/<year>/DataFiles/`
- Linked mortality: NCHS 2019 public-use LMF (follow-up through 2019-12-31),
  data dictionary dated 2026-01-14

## Commands and scripts

| Script | Command | Output |
|---|---|---|
| `scripts/01_download.py` | `python (configured project environment) analysis\01-nhanes\scripts\01_download.py` | `data/raw/nhanes/`, `manifest.csv` |
| `scripts/02_harmonize.py` | `python (configured project environment) analysis\01-nhanes\scripts\02_harmonize.py` | `data/processed/nhanes/analysis_dataset.csv` |
| `scripts/03_survey_analysis.R` | `Rscript analysis\01-nhanes\scripts\03_survey_analysis.R` | `results/tables/01-nhanes/table0鈥able7` |
| `scripts/04_rcs.R` | `Rscript analysis\01-nhanes\scripts\04_rcs.R` | `table8_rcs_cvd_by_t_score.csv` |
| `scripts/05_km.R` | `Rscript analysis\01-nhanes\scripts\05_km.R` | `table9_km_allcause_by_bmd.csv` |
| `scripts/06_figure1.py` | `python 鈥?06_figure1.py` (project `.venv`) | `figures/fig1_nhanes_bmd_cvd.png` |

`survey::svykm()` with design-based SEs did not terminate on this sample and was
replaced by the explicit weighted estimator in `05_km.R`.

## Outputs

- `data/processed/nhanes/analysis_dataset.csv` 鈥?one row per participant (14,260), 52 columns
- `data/processed/nhanes/variable_dictionary.csv` 鈥?every derived variable and its source item
- `data/processed/nhanes/t_score_reference.csv` 鈥?the young-adult reference actually used
- `results/tables/01-nhanes/table0_sample_flow.csv` 鈥?sample derivation
- `鈥?table1_weighted_by_bmd_category.csv` 鈥?weighted characteristics by BMD category
- `鈥?table2_cvd_prevalence_by_bmd.csv` 鈥?weighted CVD prevalence, overall and by sex
- `鈥?table3_survey_logistic_cvd.csv` 鈥?ORs for prevalent CVD (M1鈥揗4, subgroups, continuous, trend)
- `鈥?table4_survey_cox_mortality.csv` 鈥?HRs for all-cause and CVD mortality
- `鈥?table5_sensitivity.csv` 鈥?self-reported osteoporosis, sex-neutral T-score, lumbar spine BMD
- `鈥?table6_event_counts.csv` 鈥?event counts
- `鈥?table7_interaction_tests.csv` 鈥?BMD 脳 sex and BMD 脳 age interaction
- `鈥?table8_rcs_cvd_by_t_score.csv` 鈥?spline curve with 95% band
- `鈥?table9_km_allcause_by_bmd.csv` 鈥?weighted KM by BMD category
- `figures/fig1_nhanes_bmd_cvd.png` 鈥?four-panel figure

## Key results (as computed by these scripts)

- Weighted prevalent CVD: normal BMD 9.8%, osteopenia 13.6%, osteoporosis 22.9%;
  the crude gradient is present in both sexes.
- The adjusted association with *prevalent* CVD is weak and not robust: fully
  adjusted OR for osteoporosis vs normal 1.22 (95% CI 0.85鈥?.75), p for trend
  0.30. It is concentrated in the 40鈥?4 age band (OR 2.19, 1.26鈥?.81) and in
  stroke (OR 1.62, 1.09鈥?.39; trend p = 0.03) rather than in coronary disease
  (OR 1.11, 0.74鈥?.68). BMD 脳 sex interaction p = 0.047, BMD 脳 age p = 0.071.
- The mortality gradient is steep and survives full adjustment: osteoporosis vs
  normal HR 1.76 (1.42鈥?.18) for all-cause death, 1.51 (1.05鈥?.17) for CVD death,
  1.66 (1.23鈥?.25) in women and 1.92 (1.44鈥?.56) in men.
- Spline: monotonic gradient in prevalent CVD with no threshold (spline term
  p = 0.005), but the effect size is modest once adjusted.
- Interpretation: in NHANES the bone鈥搗ascular link behaves much more like a
  shared *prognostic* process than a shared *prevalent* disease state. This
  asymmetry is exactly what the genetic layer (module `02-mr`) has to arbitrate.

## Limitations

- Prevalent CVD is self-reported physician diagnosis, not adjudicated events;
  self-report misclassification is likely non-differential but attenuates
  estimates.
- NHANES is cross-sectional for the BMD鈥揅VD association; only the mortality
  analyses are prospective.
- Femoral neck DXA is unavailable in cycles 2011鈥?2 and 2015鈥?6, so those waves
  contribute no BMD and the pooled weights use five cycles (WTMEC2YR/5).
- Cerebrovascular deaths cannot be separated from "all other causes" in cycles
  2015鈥?8 because the public-use leading-cause recode releases only 001, 002 and
  010 for those years; `cvd_death` is therefore under-ascertained in the two most
  recent cycles.
- The young-adult reference is drawn from the same survey cycles; it matches
  published NHANES-derived references closely but is not an external standard.
- Visual inspection of `fig1` was limited to the geometric collision scan and a
  border/margin check: no vision-capable model is configured in this session.


