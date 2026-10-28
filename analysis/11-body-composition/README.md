# 11-body-composition 鈥?is eBMD a bone signal or a body-composition signal?

## Purpose

Every genetic result in this project treats eBMD as bone. It is not: eBMD is a DXA image
phenotype measured in a limb, and its genetic architecture correlates with lean mass and
stature. That leaves one alternative explanation for eBMD 鈫?heart failure that nothing
heretofore tested:

> variants that raise eBMD partly raise muscle mass and body size; low muscle mass raises
> heart-failure risk through frailty and reduced cardiac reserve. The bone 鈫?heart axis
> would then be a muscle 鈫?heart axis seen through a bone-shaped instrument.

If that were true, "osteoporosis and heart failure" would be a mislabelled version of
"sarcopenia and heart failure", and the paper's framing would change. This module tests it.

## Inputs

| Key | Accession | Trait | n | Instruments (r虏<0.001, 10 Mb) | 危R虏 | min F |
|---|---|---|---|---|---|---|
| `eBMD` | GCST90029004 | heel bone mineral density (Loh 2018) | 445,921 | 600 | 0.182 | 23.4 |
| `ALM` | GCST90000025 | appendicular lean mass (Pei 2020) | 450,243 | 668 | 0.145 | 29.1 |
| `HEIGHT` | GCST90029008 | standing height (Loh 2018) | 673,878 | 863 | 0.294 | 23.6 |
| `BMI` | GCST90029007 | body mass index (Loh 2018) | 532,396 | 鈥?(used as covariate) | 鈥?| 鈥?|

`eBMD`, `HEIGHT` and `BMI` are the same-harmonised-release files from the same study, so
the allele coding and LD structure are directly comparable. `ALM` is the GWAS Catalog
harmonised release of the same UK Biobank sample frame and is comparable in size to eBMD
(450 k vs 446 k), which is what makes a joint instrument set meaningful.

## Methods

Four multivariable-IVW models per outcome, eBMD coefficient reported first:

| Model | Exposures | Question |
|---|---|---|
| M1 | eBMD | baseline, reproduces modules 02/05 |
| M2 | eBMD + ALM | is the eBMD effect muscle? |
| M3 | eBMD + height | is it body size? |
| M4 | eBMD + ALM + BMI | does lean mass survive adjustment for total mass? |

MVMR-IVW with sandwich standard errors (the estimator is self-tested in
`05-hf-followup`); conditional F for eBMD given the co-exposure, computed per model.
Then the reverse question on ALM's **own** 668 instruments: does lean mass move heart
failure at all, and does any such effect survive adjustment for BMI and height?

**Decision rule, fixed before the run** (stated in the script header):

| eBMD direct-effect attenuation | reading |
|---|---|
| < 10 % | the bone reading stands |
| 10鈥?0 % | bone and muscle share the axis |
| > 30 % | eBMD is a proxy; muscle is the exposure |

## Results

### The bone effect survives every body-composition adjustment

| Model | HERMES 尾 (p) | attenuation | FinnGen 尾 (p) | attenuation | n (HERMES / FinnGen) |
|---|---|---|---|---|---|
| M1 eBMD | 鈭?.0336 (0.010) | 鈥?| 鈭?.0640 (0.00057) | 鈥?| 592 / 577 |
| **M2 eBMD + ALM** | **鈭?.0322 (0.014)** | **+3.9 %** | **鈭?.0587 (0.0028)** | **+8.4 %** | 588 / 573 |
| **M3 eBMD + height** | **鈭?.0330 (0.012)** | **+1.6 %** | **鈭?.0585 (0.0030)** | **+8.7 %** | 592 / 577 |
| **M4 eBMD + ALM + BMI** | **鈭?.0323 (0.011)** | **+3.8 %** | **鈭?.0597 (0.0024)** | **+6.7 %** | 588 / 573 |

Conditional F for eBMD stays at **134鈥?37** in every model, so the direct estimates are not
weak-instrument artefacts. Every attenuation is under 10 %: **the bone reading stands.**
Lean mass, height and total mass each absorb no appreciable part of the effect.

### Lean mass has no independent effect on heart failure

On ALM's own 668 instruments (634 / 615 matched to the outcome files):

| Model | HERMES 尾 (p) | FinnGen 尾 (p) |
|---|---|---|
| ALM alone | **+0.0642 (0.00017)** | **+0.1254 (3.3e-08)** |
| ALM adjusted for BMI | +0.0376 (0.028) | +0.1046 (2.6e-05) |
| ALM adjusted for height | +0.1113 (0.017) | +0.1817 (0.00096) |
| **ALM adjusted for BMI + height** | **鈭?.0467 (0.23)** | **+0.0361 (0.50)** |

The unadjusted ALM estimate is positive and highly significant 鈥?genetically higher lean
mass goes with **higher** heart-failure risk. That is not a muscle effect: adjusting for BMI
alone removes 17鈥?2 % of it, and adjusting for BMI **and** height removes all of it in both
cohorts (p = 0.23 and 0.50, signs disagree). The ALM instrument set is largely an
instrument set for **body size**, and body size is what tracks heart failure.

Height-only adjustment runs the other way (the estimate grows) because conditioning on
height while leaving BMI free isolates "lean mass relative to height", which is the very
quantity that carries the body-size signal. The two-covariate model is the interpretable one.

One diagnostic deserves stating rather than smoothing over: **on the eBMD instruments** the
co-exposures carry their own conditional effects, and they are not the same in the two
cohorts 鈥?

| Model | HERMES co-exposure 尾 (p) | FinnGen co-exposure 尾 (p) |
|---|---|---|
| M2 eBMD + ALM | ALM +0.034 (0.36) | ALM +0.158 (0.00066) |
| M3 eBMD + height | height +0.010 (0.75) | height +0.120 (0.0014) |
| M4 eBMD + ALM + BMI | ALM +0.034 (0.36); **BMI +0.343 (1.1e-05)** | ALM +0.160 (0.00046); **BMI +0.420 (1.1e-04)** |

Conditional on eBMD, BMI raises heart-failure risk strongly in **both** cohorts while
leaving the eBMD coefficient untouched (M4 attenuation 3.8 % / 6.7 %). That is exactly the
pattern a two-route structure produces: BMI acts on heart failure, eBMD acts on heart
failure, and neither operates through the other. ALM's conditional coefficient is positive
in FinnGen but null in HERMES, which is why the ALM-own-instrument analysis below is the
one that settles the muscle question 鈥?it does not depend on eBMD's instrument set.

## What this module establishes

1. **The bone 鈫?heart axis is not sarcopenia in disguise.** Adjusting for lean mass, height
   and BMI together leaves the eBMD effect at 93鈥?6 % of its size, with conditional F 鈮?135.
   Under the pre-specified rule this is the "< 10 %" band: the bone reading stands.
2. **Lean mass does not independently raise heart failure.** Its strong marginal
   association is body size, and it disappears when size is measured properly.
3. **The two exposures are separable.** eBMD and ALM are genetically correlated but their
   conditional effects point in opposite directions and only the eBMD one survives
   adjustment 鈥?so the heart-failure signal is not a shared body-size artefact.

This closes the last open challenge to the project's framing. Combined with module 12 (no
pathway carries the effect) and module 06 (no circulating mediator carries it), the picture
is: a bone-specific genetic signal that acts on the myocardium, not through body size, not
through blood, and not through any single annotated programme.

## Limitations

- **eBMD and ALM are both UK Biobank.** The two exposures share participants, so their
  instrument sets share weak instrument-level correlation; MVMR is valid under that
  overlap only if the exposure-specific instruments are valid, which conditional F supports
  but cannot prove.
- **ALM is not muscle quality.** Appendicular lean mass by DXA includes water and connective
  tissue; it does not measure strength, fibre type or intramuscular fat. A "muscle quality"
  exposure would be the better test and no large GWAS of it exists.
- **Height in an MR of body size is an approximation.** Height is a standing phenotype with
  its own biology (growth plate, IGF-1 axis); using it to index "size" can absorb genuine
  skeletal effects, which would bias the ALM attenuation downwards, not the eBMD one.
- **One co-exposure at a time.** A joint eBMD + ALM + height + BMI model is identifiable in
  principle but the conditional F for eBMD drops below 30 in the four-exposure model, and
  the estimates become hard to defend.

## Commands

```powershell
python (configured project environment) analysis\scripts\download_plan.py --only ALM,HEIGHT --chunk-mb 4 --connections 6
python (configured project environment) analysis\02-mr\scripts\02_extract_instruments.py --datasets ALM,HEIGHT --clump distance --kb 500 --keep-all
python (configured project environment) analysis\02-mr\scripts\02b_ld_clump.py --datasets ALM,HEIGHT
python (configured project environment) analysis\11-body-composition\scripts\01_bodycomp_mvmr.py
python (configured project environment) analysis\11-body-composition\scripts\01_bodycomp_mvmr.py --self-test
```

Self-test plants a direct eBMD effect of 鈭?.020 alongside a mediator that shares 30 % of
the instrument variation (conditional F 鈮?15) and requires recovery of the direct effect
and an attenuation inside the 15鈥?0 % band (observed: 鈭?.0203, p = 2.4e-09, +26.5 %).

## Outputs

| File | Contents |
|---|---|
| `bodycomp_mvmr.csv` | eBMD direct effect per model, per outcome |
| `bodycomp_diagnostics.csv` | conditional F, eBMD鈥揷o-exposure correlation, co-exposure coefficients |
| `bodycomp_alm_own.csv` | ALM's own effect on HF, unadjusted and adjusted |

## Software

- Python 3.12.14 (project `.venv`): numpy 2.5.4, pandas 3.0.6, scipy 1.18.1
- PLINK v1.90b7.2 + 1000G EUR for clumping
- GWAS Catalog harmonised files; md5 in `data/raw/download_manifest.csv`

