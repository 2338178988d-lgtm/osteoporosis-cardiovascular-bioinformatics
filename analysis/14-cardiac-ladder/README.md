# 14-cardiac-ladder 鈥?where on the path to heart failure the effect sits

## Purpose

Modules 02鈥?3 established that genetically predicted heel-ultrasound bone density lowers the
risk of heart failure, and ruled out coronary disease, circulating mediators, body size and
every annotated pathway. What they never did was **measure** where the effect acts: "the
myocardium" was reached by elimination, because FinnGen R12 publishes no "heart failure
without coronary disease" endpoint.

That endpoint does not exist, but eight others do, and each occupies a defined rung between
vascular load and myocardial failure. Running the same 600 instruments against all of them
turns an inference into a measurement 鈥?or shows that the effect cannot be localised.

## Inputs

| Path | Content | Source |
|---|---|---|
| `results/tables/02-mr/instruments_eBMD.csv` | 600 LD-independent eBMD instruments | module `02-mr` |
| `data/raw/finngen/finngen_R12_I9_*.gz.tbi` | 8 tabix indexes, 1.7 MB each | FinnGen R12, `download_plan.py --category ladder` |
| `data/processed/mr/outcome_snps/FG_*.csv` | instrument-SNP rows, cached after first fetch | written by `01_ladder.py` |
| `results/tables/02-mr/instruments_eBMD.csv` 鈫?`data/processed/mr/outcome_snps/FG_HF.csv` | heart-failure anchor | module `02-mr` / `05-hf-followup` |

The 815 MB summary-statistics file for each endpoint is **never downloaded.** FinnGen
publishes a tabix index beside every file, so `analysis/scripts/_tabix_http.py` reads the
~600 instrument positions over HTTP range requests: about 30 MB per endpoint instead of
815 MB. Fetched rows are cached and the requested variant list is snapshotted in a `.snps`
file, so a rerun is offline and does not re-fetch.

## Methods

- **Exposure.** UK Biobank heel-ultrasound eBMD, GCST90029004 (Morris 2018), 600 LD-independent
  instruments (PLINK r虏 < 0.001 within 10 Mb, 1000G EUR), 危R虏 = 0.182, min F = 23.4. Identical
  to modules 02鈥?3; no instrument is re-selected per outcome.
- **Outcomes.** Eight FinnGen R12 endpoints, all Finnish, so there is no UK Biobank overlap
  with the exposure:

  | Key | Phenocode | Endpoint | Cases |
  |---|---|---|---|
  | `FG_HYPERTENS` | `I9_HYPTENS` | hypertension | 154,630 |
  | `FG_MI_STRICT` | `I9_MI_STRICT` | myocardial infarction, strict | 31,666 |
  | `FG_AF` | `I9_AF` | atrial fibrillation and flutter | 63,532 |
  | `FG_VHD` | `I9_VHD` | valvular heart disease | 28,602 |
  | `FG_CAVS` | `I9_CAVS_OPERATED` | calcific aortic stenosis, operated | 12,418 |
  | `FG_HYPERTENS_HD` | `I9_HYPERTENSHD` | hypertensive heart disease | 11,728 |
  | `FG_NONISCH_CM` | `I9_NONISCHCARDMYOP` | non-ischaemic cardiomyopathy | 13,917 |
  | `FG_CARDMYO` | `I9_CARDMYO` | cardiomyopathy, all | 7,649 |
  | `FG_HF` (anchor) | `I9_HEARTFAIL` | heart failure, strict | 37,653 |

- **Estimator.** IVW with random effects (Cochran's Q from the same fit), MR-Egger intercept as
  a directional-pleiotropy check, and leave-one-out over all instruments. Allele harmonisation
  drops palindromic variants with intermediate frequency. Unchanged from module `13-robustness`.
- **Multiplicity.** Eight endpoints are tested (the heart-failure anchor is reported in the
  same table but not counted): Bonferroni 伪 = 0.05 / 8 = 0.00625, reported alongside the
  unadjusted p.
- **Power.** Reported per endpoint as the power to detect |尾| = 0.064, the FinnGen heart-failure
  estimate, at 伪 = 0.05; and the minimum detectable effect at 80 % power.
- **Decision rules were fixed in the script header before the run** and are printed by
  `verdict()`. `02_ladder_concordance.py` then asks whether the protective rungs are one
  signal or several, by correlating per-variant Wald ratios between endpoints.

## Results

| Endpoint | 尾 | 95 % CI | p | power | Egger p | Bonferroni |
|---|---|---|---|---|---|---|
| Hypertension | +0.008 | 鈭?.032, +0.048 | 0.71 | 88 % | 0.49 | 鈥?|
| Myocardial infarction | +0.007 | 鈭?.043, +0.057 | 0.79 | 71 % | 0.06 | 鈥?|
| **Atrial fibrillation** | **鈭?.086** | 鈭?.138, 鈭?.033 | **0.0013** | 67 % | 0.73 | **yes** |
| **Valvular heart disease** | **鈭?.074** | 鈭?.119, 鈭?.030 | **0.0012** | 80 % | 0.72 | **yes** |
| Calcific aortic stenosis | 鈭?.024 | 鈭?.090, +0.042 | 0.47 | 47 % | 0.11 | 鈥?|
| Hypertensive heart disease | 鈭?.013 | 鈭?.083, +0.056 | 0.70 | 44 % | 9e-05 | 鈥?|
| **Non-ischaemic cardiomyopathy** | **鈭?.101** | 鈭?.153, 鈭?.050 | **1.0e-04** | 69 % | 0.07 | **yes** |
| Cardiomyopathy (all) | 鈭?.062 | 鈭?.132, +0.008 | 0.085 | 43 % | 0.009 | 鈥?|
| Heart failure (anchor) | 鈭?.064 | 鈭?.100, 鈭?.028 | 0.00057 | 93 % | 0.80 | 鈥?|

Leave-one-out: no significant estimate crosses zero for any single-variant deletion
(heart failure 鈭?.068 to 鈭?.059; atrial fibrillation 鈭?.094 to 鈭?.080; valve disease 鈭?.081 to
鈭?.069; non-ischaemic cardiomyopathy 鈭?.109 to 鈭?.089).

### What the pre-registered rules said

- **Rule A** (non-ischaemic cardiomyopathy *and* hypertensive heart disease protective) 鈥?
  **not met**: non-ischaemic cardiomyopathy is strongly protective, hypertensive heart disease
  is null at 44 % power.
- **Rule B** (only heart failure survives) 鈥?**not met**: three rungs survive Bonferroni.
- **Rule C** (calcific aortic stenosis protective) 鈥?**not met**: 鈭?.024, p = 0.47, 47 % power.

None of the three pre-registered outcomes occurred. That is reported as such; the rules are
not reinterpreted after the fact.

### The pattern that was observed

Every **non-coronary** cardiac endpoint has a negative point estimate; the three with adequate
power are Bonferroni-significant. The two rungs that are truly null and adequately powered are
**hypertension** (88 % power) and **myocardial infarction** (71 % power) 鈥?the two that would
have to move if the effect acted through blood pressure or atherosclerosis. The three
non-coronary nulls (calcific aortic stenosis, hypertensive heart disease, cardiomyopathy-all)
are the three lowest-powered rungs (43鈥?7 %) and all point the same way.

So the effect is not "heart-failure specific", as the previous version of this project
concluded; it is **heart-failure specific only because heart failure was the only non-coronary
endpoint anyone had tested.**

### The three numbers are one number

If heart failure is a mixture of an ischaemic part (where the effect is absent) and a
non-ischaemic part (where it is present), the all-HF and HF+CHD estimates together imply what
the non-ischaemic part must be. `03_dilution.py` does that arithmetic with the mixing weight
w = 17,532 / 37,653 = 0.466:

| Quantity | b (95 % CI) | p |
|---|---|---|
| all heart failure (measured) | -0.064 (-0.100, -0.028) | 0.00057 |
| heart failure + coronary disease (measured) | +0.008 (-0.046, +0.062) | 0.78 |
| **heart failure without coronary disease (implied)** | **-0.127 (-0.209, -0.044)** | 0.0028 |
| non-ischaemic cardiomyopathy (measured separately) | -0.102 (-0.153, -0.050) | 1.0e-04 |
| implied minus observed | +0.025 (se 0.050) | **0.61** |

The separately measured non-ischaemic cardiomyopathy estimate and the value implied by the
mixture are statistically indistinguishable. **The ladder is one effect diluted by the
ischaemic half of heart failure**, not three findings. That form is falsifiable: any new
non-coronary cardiac endpoint should land between -0.10 and -0.13, and any aggregate that
mixes in ischaemic cases should be diluted to -0.05 to -0.07.

The arithmetic assumes the two subtypes partition heart failure and that the effects are
additive on the same scale. FinnGen constructs its subtypes against different control sets, so
they do not partition exactly; the calculation is an approximation and is labelled as one.

### One signal, observed several times

The Finnish endpoints are ICD-code sets that share patients, so three concordant estimates are
not three independent confirmations. `02_ladder_concordance.py` correlates per-variant Wald
ratios across endpoints:

|  | HF | AF | VHD | NONISCH_CM | HTN | MI |
|---|---|---|---|---|---|---|
| HF | 1.00 | 0.53 | 0.41 | 0.69 | 0.29 | 0.29 |
| AF | 0.53 | 1.00 | 0.55 | 0.42 | 0.50 | 0.31 |
| VHD | 0.41 | 0.55 | 1.00 | 0.40 | 0.29 | 0.15 |
| NONISCH_CM | 0.69 | 0.42 | 0.40 | 1.00 | 0.11 | 0.02 |
| HTN | 0.29 | 0.50 | 0.29 | 0.11 | 1.00 | 0.31 |
| MI | 0.29 | 0.31 | 0.15 | 0.02 | 0.31 | 1.00 |

Mean r within the protective rungs 0.50; between protective and null rungs 0.25. The same
variants carry the effect across the protective endpoints 鈥?one axis, read several times 鈥?
and the weakest agreement in the table is non-ischaemic cardiomyopathy against myocardial
infarction (r = 0.02).

## D. The discriminating panel: three questions the ladder left open

Eight endpoints located the effect but did not say what kind of thing it is. Twelve further
endpoints, all read through their tabix indexes, settle three of those questions
(`04_discriminating_panel.py`). Realised power is lower than a case-count projection suggests,
because Cochran's Q is significant for these endpoints and the random-effects standard error
inflates accordingly; the power column below is the realised value.

### D1. Cardiac, not systemic matrix

The competing reading after the first eight endpoints was that the exposure tracks
collagen/matrix quality everywhere, and that heart, veins and aorta should all move. They do not.

| Endpoint | b | 95 % CI | p | power |
|---|---|---|---|---|
| varicose veins (38,467) | -0.014 | -0.068, +0.041 | 0.63 | 63 % |
| venous thromboembolism (26,333) | -0.049 | -0.100, +0.002 | 0.061 | 69 % |
| pulmonary embolism (12,762) | -0.038 | -0.098, +0.023 | 0.22 | 54 % |
| aortic aneurysm (9,937) | -0.027 | -0.095, +0.041 | 0.44 | 45 % |

None survives the family threshold (0.0125). **Rule A fires: the effect is cardiac-specific and
the systemic-matrix reading is excluded.** One caveat is recorded rather than smoothed over:
venous thromboembolism sits at -0.049, about half the cardiac effect, and after removing four
heart-failure-adjacent loci it reaches p = 0.046. It is the one row that is not a clean null,
and it is a thrombosis endpoint where the pathway need not run through matrix at all.

### D2. Atrial fibrillation, not the conduction system, and not arrhythmia in general

| Endpoint | b | 95 % CI | p | power |
|---|---|---|---|---|
| **other arrhythmias (39,265)** | -0.026 | -0.060, +0.009 | 0.14 | **95 %** |
| paroxysmal tachycardia (12,878) | -0.031 | -0.088, +0.027 | 0.29 | 59 % |
| conduction disorders (14,027) | -0.019 | -0.071, +0.033 | 0.47 | 67 % |
| AV block (7,850) | -0.061 | -0.129, +0.008 | 0.084 | 45 % |
| **atrial fibrillation, reimbursement (40,251)** | **-0.096** | -0.155, -0.037 | **0.0015** | 56 % |
| atrial fibrillation, strict (63,532) | -0.085 | -0.138, -0.033 | 0.0013 | 67 % |

**Rule B fires: atrial moves, conduction does not.** Two points carry the weight. "Other
arrhythmias" is the non-AF arrhythmia bucket at **95 % power** and is null, so the effect is not
generic arrhythmia. Conduction disorders are null at 67 % power, so it is not the fibrous
skeleton of the heart either. The reimbursement definition reproduces at -0.096 - a different
ascertainment of the same arrhythmia, and *larger* than the strict one.

### D3. Aggregates carry no effect

| Endpoint | b | p | power |
|---|---|---|---|
| cardiovascular disease, all (247,538) | -0.016 | 0.30 | 99 % |
| ischaemic heart disease (84,088) | +0.015 | 0.42 | 94 % |
| major CHD event (56,650) | +0.029 | 0.16 | 88 % |
| heart failure (37,653) | -0.064 | 0.00057 | 93 % |
| non-ischaemic cardiomyopathy (13,917) | -0.102 | 1.0e-04 | 69 % |

Strict monotonicity **fails**: the three aggregates sit within one standard error of zero and
their signs are noise (the first two differ by 0.0008 against a standard error of about 0.02).
The informative test is the pooled contrast, which the figure draws:

- aggregates pooled (all CVD + IHD + major CHD): **+0.005 (se 0.010)**
- specific pooled (heart failure + non-ischaemic cardiomyopathy): **-0.077 (se 0.015)**
- difference **-0.081 (se 0.018), p = 7.9e-06**

**Rule C fires in its usable form.** The effect appears below the level at which endpoints
aggregate, which is what the mixture synthesis in "The three numbers are one number" requires.
The strict monotonicity test is reported as failed, not discarded.

### D4. Whether the valvular effect is calcific stays open

FinnGen publishes only three valve endpoints and no mitral or aortic-regurgitation codes, so
"valvular disease is protective" could not be split by lesion type from Finnish data. A direct
measurement of aortic valve calcification (GCST90310293, 14,819 cases / 941,863 controls,
PMID 38494474) was added (`06_valve_calcification.py`).

| Endpoint | b | 95 % CI | p | power |
|---|---|---|---|---|
| valvular heart disease, ICD aggregate | -0.074 | -0.119, -0.030 | 0.0012 | 80 % |
| calcific aortic stenosis, operated | -0.024 | -0.090, +0.042 | 0.47 | 47 % |
| **aortic valve calcification** | **-0.058** | -0.119, +0.004 | 0.065 | 53 % |

Direct valve calcification is directionally consistent with the aggregate and about
three-quarters of its size, but it does not reach significance and its CI spans the aggregate
estimate. **Rule D returns "open", not "negative"** - the question needs a larger calcification
GWAS, and the honest statement is that the valvular effect is compatible with a calcific
component without establishing one.

## E. Checks a reviewer will ask for

`05_panel_checks.py` recomputes every estimate three ways. None of the panel conclusions moves.

| Check | Result |
|---|---|
| **Shared-locus removal** - drop eBMD instruments within 500 kb of any of the 754 genome-wide-significant heart-failure variants | 4 loci removed; every estimate essentially unchanged (AF -0.086 to -0.094; non-ischaemic CM -0.102 to -0.103; VHD -0.075 to -0.076). One estimate changes sign (MI, +0.007 to -0.001) and it is a null either way |
| **Steiger filtering** - drop variants explaining more variance in the outcome than the exposure, using module 05's estimator | 0-3 variants removed per endpoint; no estimate moves by more than 0.003 |
| **Effective number of tests** - Li & Ji (2005), from the eigenvalue spectrum of the 20-endpoint Wald-ratio correlation matrix | Meff = **14.1** for 20 endpoints, so the corrected threshold is 0.0035 rather than Bonferroni's 0.0025. AF (0.0013), AF-reimbursement (0.0015), VHD (0.0012) and non-ischaemic CM (1.0e-04) pass under either |

The shared-locus check exists partly because an earlier version of it silently did nothing: it
queried a tabix index for `I9_HEARTFAIL` that was never downloaded (only the full `.gz` was), so
every "-N loci" in the first run was actually zero. It now reads the local file and reports 754
top hits. The bug is described in the function docstring so it is not reintroduced.

## G. Five more outcome cohorts: does it replicate outside Finland?

Sections A-F are Finnish except the heart-failure anchor. Two objections follow: atrial
fibrillation - the largest and most surprising non-coronary effect - was measured once, in one
founder population; and the heart-failure anchor that did replicate shares UK Biobank
participants with the exposure. Five further outcome cohorts settle both
(`07_big_cohorts.py`).

Four of them publish a tabix index beside the harmonised file, so the 600 instrument positions
are read over HTTP range requests and the 1.7-4.0 GB files are never downloaded. The fifth
(GCST90475970) has no index and was fetched whole, 579 MB.

| key | source | cases | ancestry | overlaps |
|---|---|---|---|---|
| `AF_US` | GCST90475970, PheCode 427.21, MVP (PMID 39024449) | 71,577 | US European | **neither UKB nor FinnGen** |
| `AF_NFE` | GCST90473571, ICD10 I48, UKB (PMID 40770095) | 40,273 | non-Finnish European | UKB only |
| `AF_BIG` | GCST90624411 (PMID 40645996) | 252,438 | European | UKB + FinnGen |
| `HF_BIG` | GCST90455657 (PMID 40195560) | 168,363 | European | UKB + FinnGen |
| `HF_NFE` | GCST90473575, ICD10 I50, UKB (PMID 40770095) | 19,629 | non-Finnish European | UKB only |

### Atrial fibrillation

| cohort | b | 95 % CI | p | power |
|---|---|---|---|---|
| FinnGen reimbursement (40,251) | -0.096 | -0.155, -0.037 | 0.0015 | 56 % |
| FinnGen strict (63,532) | -0.086 | -0.138, -0.033 | 0.0013 | 67 % |
| **GCST90475970, US (71,577)** | **-0.071** | **-0.116, -0.025** | **0.0023** | 79 % |
| GCST90624411 (252,438) | -0.055 | -0.091, -0.020 | 0.0021 | 95 % |
| GCST90473571, non-Finnish (40,273) | -0.033 | -0.076, +0.011 | 0.14 | 82 % |
| **pooled, k=5** | **-0.063** | **-0.083, -0.042** | **1.7e-09** | 鈥?|

Five cohorts, five negative estimates, I^2 = **0.03**. **Rule 1 fires: the effect reproduces in
a cohort that overlaps neither UK Biobank nor FinnGen** (-0.071, p = 0.0023). Atrial
fibrillation is no longer a Finnish-only result.

### Heart failure

| cohort | b | 95 % CI | p | power |
|---|---|---|---|---|
| FinnGen strict (37,653) | -0.064 | -0.100, -0.028 | 0.00057 | 93 % |
| GCST90473575, non-Finnish (19,629) | -0.051 | -0.097, -0.004 | 0.032 | 78 % |
| HERMES (47,309) | -0.034 | -0.059, -0.008 | 0.010 | 100 % |
| GCST90455657 (168,363) | -0.020 | -0.044, +0.004 | 0.11 | 100 % |
| **pooled, k=4** | **-0.037** | **-0.056, -0.019** | **1.0e-04** | 鈥?|

Four cohorts, four negative estimates, I^2 = 0.32. No cohort reverses. The 168,363-case file is
the *weakest* of the four (-0.020) at 100 % power, which is the opposite of what winner's curse
would produce in a single study and is consistent with a real but small effect.

### The Finnish estimates are the largest in both panels

- atrial fibrillation: Finnish pooled -0.090 (k=2) vs other -0.053 (k=3), difference -0.037,
  p = 0.11
- heart failure: Finnish -0.064 (k=1) vs other -0.029 (k=3), difference -0.035, p = 0.09

Neither difference reaches significance, and both point the same way. With five and four
cohorts that is not enough to establish a founder-population effect, but it is enough to say
the *pooled* estimate should be quoted rather than the FinnGen one. The conservative reading is
that the true effect is nearer -0.05 for atrial fibrillation and -0.035 for heart failure than
the -0.086 and -0.064 that FinnGen alone gives.

### What this adds, and what it does not

- **Adds**: the atrial-fibrillation effect and the heart-failure effect both survive cohorts
  that use different case definitions (PheCode, ICD10, FinnGen strict), different ancestries
  within Europe, and different degrees of participant overlap with the exposure. No cohort
  reverses sign.
- **Does not add**: an East Asian or African replication (module 08's BioBank Japan null is
  untouched and remains underpowered), or a cohort independent of UK Biobank for heart failure
  (none of the four HF cohorts is free of UKB *and* FinnGen; only atrial fibrillation has that).
- **A technical note**: GCST90475970 publishes an odds ratio with `standard_error` set to `#NA`
  and a confidence interval instead. It is converted to a log-OR beta with
  SE = (ln(ci_upper) - ln(ci_lower)) / 3.92, which reproduces the stated OR and interval
  exactly. Without that conversion the file yields no usable rows at all - the first run of
  this script silently produced nothing for that cohort.

## H. Non-ischaemic cardiomyopathy does not replicate, and that is the finding

Section G replicated atrial fibrillation outside Finland. The largest effect in the whole
project - non-ischaemic cardiomyopathy, -0.102 - had never been tested outside Finland. Three
cardiomyopathy files exist with public summary statistics, and all three were run
(`07_big_cohorts.py`).

| cohort | b | 95 % CI | p | power | case definition |
|---|---|---|---|---|---|
| FinnGen non-ischaemic (13,917) | **-0.102** | -0.153, -0.050 | 1.0e-04 | 69 % | ICD non-ischaemic CM |
| GCST90473562, non-Finnish (3,578) | -0.070 | -0.160, +0.021 | 0.13 | **28 %** | ICD cardiomyopathy |
| GCST90475951, MVP (19,708) | -0.001 | -0.053, +0.052 | 0.98 | 67 % | PheCode cardiomyopathy |
| GCST90475949, MVP (20,680) | +0.003 | -0.048, +0.054 | 0.90 | 69 % | PheCode cardiomyopathy |
| **pooled non-Finnish, k=3** | **-0.009** | -0.043, +0.025 | 0.62 | 鈥?| I^2 = 0.01 |

**Pre-registered rule fires: both well-powered files are null (67 % and 69 % power, estimates
+0.003 and -0.001), so the cardiomyopathy effect is demoted from main evidence to a Finnish
single-cohort result.** The Finnish estimate differs from the non-Finnish pooled estimate by
-0.093 (p = 0.003).

### The caveat that must travel with this

The case definitions are **not the same**, and this is the main alternative explanation:

- FinnGen `I9_NONISCHCARDMYOP` is a narrow ICD category ("non-ischaemic cardiomyopathy").
- The two MVP files are PheCode "cardiomyopathy", which is broader and includes hypertrophic
  and alcoholic cardiomyopathy.

A broader definition should attenuate a subtype-specific effect, but attenuation and a point
estimate of +0.003 are different things; dilution alone does not usually produce an exact zero.
The one file with a comparable ICD basis (GCST90473562) gives -0.070 in the same direction and
is the weakest test here at 28 % power. So the honest reading is:

- **the narrow Finnish estimate is not confirmed anywhere else**, and
- **it is not excluded either**, because no well-powered file uses the same narrow definition.

The paper should not claim a replicated cardiomyopathy effect, and should not claim it is
refuted.

### The pattern behind it: the Finnish estimates are larger everywhere

Pooling the cohorts on one axis shows something the endpoint-by-endpoint tables hide. For every
endpoint measured in both settings, the Finnish estimate is the larger one:

| endpoint | Finnish | k | other | k | ratio | difference | p |
|---|---|---|---|---|---|---|---|
| atrial fibrillation | -0.090 | 2 | -0.053 | 3 | 1.7 | -0.037 | 0.11 |
| heart failure | -0.064 | 1 | -0.029 | 3 | 2.2 | -0.035 | 0.09 |
| non-ischaemic cardiomyopathy | -0.102 | 1 | -0.009 | 3 | **11.8** | -0.093 | **0.003** |

Three of three point the same way; one reaches significance.

`08_ascertainment_check.py` rules out one of the three candidate explanations. If Finnish case
ascertainment inflated *everything* - a pipeline or coding effect - then Finnish endpoints that
should not respond to bone density would be inflated too. They are not:

| Finnish estimates | mean abs(b) | largest / smallest abs(b) | median abs(z) |
|---|---|---|---|
| six negative controls | 0.015 | **largest 0.034** | 0.55 |
| four protective endpoints | 0.087 | **smallest 0.064** | 3.33 |

The two groups do not overlap, by a factor of 1.9. A general ascertainment or pipeline inflation
is therefore not visible at this precision. What remains is either an inflation that acts only on
endpoints with a real effect (registry coding admits severe cases and omits mild ones, which
inflates a per-SD estimate) or a founder-population difference. Public data cannot separate those
two.

A weak additional signal points at severity: within FinnGen, atrial fibrillation defined by
reimbursement (40,251 cases, which requires treatment or admission) gives -0.096 and the strict
definition (63,532 cases) gives -0.086. The heavier definition is larger, in the direction
severity-selection predicts, but by only 0.010 - far inside one standard error.

Either way the rule for the paper holds: **quote the pooled estimate with the Finnish cohorts
marked, never the FinnGen number alone.**

One correction worth recording: the first version of this check counted the heart-failure
positive control in module 13's table as a negative control. That put the very effect under study
on both sides of the comparison and made the null bound come out at exactly 0.064, the
heart-failure estimate itself - which read as "general inflation cannot be excluded". Excluding
it changes the bound to 0.034 and reverses that conclusion.

### Consequence for what gets written

- Atrial fibrillation: 5 cohorts, pooled -0.063, reproduces in a cohort overlapping neither UK
  Biobank nor FinnGen. **Main evidence.**
- Heart failure: 4 cohorts, pooled -0.037, all negative. **Main evidence.**
- Non-ischaemic cardiomyopathy: 1 Finnish cohort, not confirmed and not excluded. **Limitation,
  with the case-definition mismatch stated.**
- Valve disease: 1 Finnish cohort, no replication data exists at all (section D4). **Limitation.**

## F. The observed power is lower than the projection

The design stage projected powers from case counts alone (varicose veins 0.94, VTE 0.82,
AF-reimbursement 0.94). The realised values are 0.63, 0.69 and 0.56. The reason is that
Cochran's Q is significant for these endpoints and the random-effects standard error inflates
by sqrt(Q/df), so a case-count projection overstates power whenever the instrument set has
heterogeneous effects - which this one does, by design (module 14's whole finding is that the
effect is dispersed). Every power figure quoted in this README is the realised one from the
fitted model, not the projection.

## Software and data sources

- Python 3.12.14 (project `.venv`): numpy 2.5.4, pandas 3.0.6, scipy 1.18.1, matplotlib 3.11.2
- FinnGen R12 release (`finngen-public-data-r12`, accessed 2026-10-23); case/control counts from
  `data/processed/finngen_R12_manifest.tsv`
- No PLINK or external LD call in this module: the instrument set is reused from module 02

## Commands and scripts

```powershell
# tabix indexes only (1.7 MB per endpoint) 鈥?the 815 MB files are not needed
python (configured project environment) analysis\scripts\download_plan.py --category ladder

# self-test, then the ladder, then the concordance check
python (configured project environment) analysis\14-cardiac-ladder\scripts\01_ladder.py --self-test
python (configured project environment) analysis\14-cardiac-ladder\scripts\01_ladder.py --threads 10
python (configured project environment) analysis\14-cardiac-ladder\scripts\02_ladder_concordance.py
python (configured project environment) analysis\14-cardiac-ladder\scripts\03_dilution.py

# discriminating panel (12 endpoints, tabix indexes only) + the three checks
python (configured project environment) analysis\14-cardiac-ladder\scripts\04_discriminating_panel.py --self-test
python (configured project environment) analysis\14-cardiac-ladder\scripts\04_discriminating_panel.py --threads 12
python (configured project environment) analysis\14-cardiac-ladder\scripts\05_panel_checks.py --self-test
python (configured project environment) analysis\14-cardiac-ladder\scripts\05_panel_checks.py

# valve: direct calcification endpoint (GCST90310293)
python (configured project environment) analysis\scripts\download_plan.py --only AVC --chunk-mb 4 --connections 8
python (configured project environment) analysis\14-cardiac-ladder\scripts\06_valve_calcification.py --self-test
python (configured project environment) analysis\14-cardiac-ladder\scripts\06_valve_calcification.py

# five further outcome cohorts (four tabix-only, one full download)
python (configured project environment) analysis\scripts\add_plan_rows_big.py
python (configured project environment) analysis\scripts\download_plan.py --only AF_US --chunk-mb 4 --connections 8

# cardiomyopathy replication inputs (section H): two tabix indexes + one 586 MB full file
python (configured project environment) analysis\scripts\add_plan_rows_cm.py
python (configured project environment) analysis\scripts\download_plan.py --only CM_MVP1,CM_MVP2,CM_NFE --chunk-mb 4 --connections 8
python (configured project environment) analysis\14-cardiac-ladder\scripts\07_big_cohorts.py --self-test
python (configured project environment) analysis\14-cardiac-ladder\scripts\07_big_cohorts.py --threads 24
python (configured project environment) analysis\14-cardiac-ladder\scripts\08_ascertainment_check.py

# figures
python (configured project environment) analysis\scripts\make_figure9.py
python (configured project environment) analysis\scripts\make_figure11.py
python (configured project environment) analysis\scripts\make_figure12.py
```

`--only FG_AF,FG_VHD` runs a subset. `--force-fetch` ignores the cache and re-reads the remote
file. Self-test plants a 鈭?.064 effect in a 592-variant set and checks it is recovered while a
null stays null (recovered 鈭?.052, p = 9.3e-05; null p = 0.29).

## Outputs

| File | Contents |
|---|---|
| `results/tables/14-cardiac-ladder/ladder.csv` | one row per endpoint: 尾, se, CI, p, Q, Egger intercept p, power, MDE, Bonferroni verdict |
| `results/tables/14-cardiac-ladder/ladder_loo.csv` | leave-one-out range per endpoint |
| `results/tables/14-cardiac-ladder/rung_concordance.csv` | pairwise per-variant Wald-ratio correlations |
| `results/tables/14-cardiac-ladder/dilution.csv` | the mixture synthesis: implied vs measured non-ischaemic heart failure |
| `results/tables/14-cardiac-ladder/panel.csv` | the twelve discriminating endpoints, with group, power and verdict |
| `results/tables/14-cardiac-ladder/dilution_series.csv` | aggregates vs specific endpoints, ordered |
| `results/tables/14-cardiac-ladder/panel_checks.csv` | base, shared-locus-removed and Steiger-filtered estimates for all 20 endpoints |
| `results/tables/14-cardiac-ladder/effective_tests.csv`, `endpoint_correlation.csv` | Meff and the 20-endpoint Wald-ratio correlation matrix |
| `results/tables/14-cardiac-ladder/avc.csv` | aortic valve calcification next to the two Finnish valve rungs |
| `results/tables/14-cardiac-ladder/big_cohorts.csv` | the five large-cohort estimates with ancestry and overlap |
| `results/tables/14-cardiac-ladder/big_cohorts_pooled.csv` | within-endpoint pooling of the new cohorts only |
| `results/tables/14-cardiac-ladder/big_cohorts_all.csv` | every cohort on one axis, Finnish vs other |
| `results/tables/14-cardiac-ladder/finnish_inflation.csv` | Finnish vs non-Finnish pooled estimate, all three endpoints |
| `results/tables/14-cardiac-ladder/ascertainment_check.csv` | Finnish negative controls vs Finnish protective endpoints, as an upper-bound check on coding inflation |
| `figures/fig9_cardiac_ladder.png` | the ladder; informative and uninformative nulls drawn differently |
| `figures/fig11_discriminating_panel.png` | the four decisions: matrix, rhythm, aggregation, valve |
| `figures/fig12_cohort_replication.png` | atrial fibrillation across five cohorts, heart failure across four |

## Limitations

- **The endpoints are ICD-code sets, not physiological categories.** Atrial fibrillation, valve
  disease, cardiomyopathy and heart failure share patients. The concordance analysis shows they
  share per-variant effects, which is the honest reading: one axis, several readouts. They are
  **not** four independent replications.
- **Hypertensive heart disease and calcific aortic stenosis are the two non-coronary rungs that
  did not move**, and they are also the two lowest-powered (44 % and 47 %). Their nulls
  constrain little; the estimates are compatible with the heart-failure effect.
- **`I9_HYPTENSHD` has an MR-Egger intercept p = 9e-05** and `I9_CARDMYO` p = 0.009. Both are
  null estimates, so the pleiotropy does not create a positive finding, but neither endpoint
  should be read as a clean null either.
- **No multivariable adjustment against coronary disease is possible here** at the level the
  subtype analysis used, because the ladder has no individual-level data; the coronary rung
  (myocardial infarction) is the control, and it is null.
- **Finnish population only.** Everything in this module is one founder population; module 08's
  East Asian null and module 13's DXA problem are unaffected by it.
- **The extended panel inherits the same limitation.** Twelve more endpoints sharpen what the
  effect is, but they are the same cohort, the same instruments and the same ICD-coding
  practice. Nothing in section D is an independent replication of section D's own findings.
- **Venous thromboembolism is the one panel row that is not a clean null** (-0.049, p = 0.061;
  p = 0.046 after shared-locus removal). It does not survive the family threshold and its point
  estimate is half the cardiac effect, but it is recorded rather than rounded to zero.
- **The valve question is not closed.** Aortic valve calcification (-0.058, p = 0.065) is
  compatible with the aggregate valve estimate and compatible with zero. Rule D returns
  "insufficient evidence", not "no calcific component".
- **Steiger filtering and shared-locus removal are near-identity checks here.** They remove at
  most four variants from 600 and change no conclusion. That is the expected result for an
  exposure whose effect is dispersed over hundreds of loci, but it also means these two checks
  have little power to detect anything subtle.
- **No well-powered cardiomyopathy file uses the same case definition as FinnGen.** The MVP
  PheCode files are broad (they include hypertrophic and alcoholic cardiomyopathy); FinnGen's
  `I9_NONISCHCARDMYOP` is narrow. The non-replication in section H is therefore a null in a
  *different* phenotype, and the narrow Finnish result stays formally untested outside Finland.
- **No cohort is independent of both UK Biobank and FinnGen for heart failure.** Only atrial
  fibrillation has such a cohort (GCST90475970). For heart failure, "replicates outside
  Finland" always means "replicates in a cohort that still shares UK Biobank participants".
- **Five and four cohorts are too few to establish a Finnish founder-population effect**, even
  though the Finnish estimates are the largest in both panels. The difference is reported with
  its p-value (0.11 and 0.09) and not interpreted further.
- **Meff = 14.1 was computed over all 20 endpoints together**, mixing the panel and the ladder.
  Per-family effective tests would be slightly smaller; the corrected thresholds would be
  slightly more generous. Reported as one pooled family, which is the conservative choice.

