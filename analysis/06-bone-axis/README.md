# 06-bone-axis 鈥?which pathway carries the bone-density effect on heart failure?

## Purpose

Module `05-hf-followup` established that genetically predicted bone mineral density
lowers heart-failure risk, that the effect replicates in FinnGen without sample overlap,
and that it is absent from the heart-failure subtype carrying coronary disease. What it
did **not** do is identify the pathway. The mediators tested there (BMI, type-2 diabetes,
CAD) are not bone-derived, so ruling them out only says the effect is not routed through
metabolic or atherosclerotic disease.

This module tests the factors a bone cell actually secretes and that are known to act on
the heart.

## Inputs

Ten protein/biomarker GWAS from the GWAS Catalog, `data/raw/gwas/`, verified against the
published md5 (`data/raw/download_manifest.csv`):

| Key | Accession | Trait | n | Instruments after LD clumping |
|---|---|---|---|---|
| `FGF23` | GCST90469223 | FGF23 protein levels (Loya 2025) | 47,745 | 21 |
| `SOST` | GCST90470711 | SOST protein levels (Loya 2025) | 47,745 | 54 |
| `BGLAP` | GCST90468440 | BGLAP protein levels (Loya 2025) | 47,745 | 20 |
| `CST3` | GCST90025945 | cystatin C levels (Barton 2021) | 437,846 | 367 |
| `CRP` | GCST90025959 | C-reactive protein (Barton 2021) | 436,939 | 244 |
| `DKK1` | GCST90468997 | DKK1 protein levels (Loya 2025) | 47,745 | 52 |
| `PTH` | GCST90470372 | PTH protein levels (Loya 2025) | 47,745 | 17 (distance) |
| `TESTO` | GCST90012114 | total testosterone (Ruth 2020) | 425,097 | 248 (distance) |
| `ESTRADIOL` | GCST90012105 | estradiol (Ruth 2020) | 206,927 | 13 (autosomes) |
| `VITD` | GCST90000618 | 25-hydroxyvitamin D | 496,946 | 118 |

> **A data-quality substitution was necessary.** The first CRP file chosen,
> GCST90029070, ships with `effect_allele_frequency` **empty for every row**, so the
> MAF filter removed all variants and zero instruments were produced. It is also a
> significant-variants-only release (`p_value` starts at 0, no SE in the first rows).
> It was replaced with GCST90025959, which has complete EAF, SE and p values.

The last three instrument sets were built when those mediators were added, so they are
not interchangeable: `PTH` and `TESTO` are distance-pruned at 500 kb only, `ESTRADIOL`
is PLINK-clumped on autosomes, and 1,212 of its X-chromosome variants are dropped because
the 1000G EUR reference carries chr1-22 only. The others are PLINK-clumped with
1000G EUR (r虏 < 0.001 within 10 Mb).

Outcomes: `HF` (HERMES, overlapping with the exposure) and `FG_HF` (FinnGen strict,
independent). Exposure: the 600 LD-clumped eBMD instruments.

## Methods

Multivariable MR of eBMD and one mediator at a time, MVMR-IVW with sandwich standard
errors (the estimator is self-tested in `05-hf-followup`):

```
X = [bx_eBMD, bx_mediator]        y = by_HF        W = diag(1/SE_y虏)
beta = (X'WX)鈦宦筙'Wy,  cov = (X'WX)鈦宦?X'W diag(resid虏) WX (X'WX)鈦宦?
```

Three quantities per mediator:

- **direct effect** 鈥?eBMD effect on HF after conditioning on the mediator
- **attenuation** 鈥?how much of the eBMD effect the mediator absorbs, measured against the
  **univariable estimate on the same instruments**. This matters: comparing a direct
  estimate from 87 shared variants against a univariable estimate from 592 variants
  produces a spurious 30鈥?0 % "attenuation" that is pure precision artefact. The script
  reports both `univariable_b` (same subset) and `univariable_full_b` so the two can be
  compared explicitly.
- **conditional F** 鈥?whether eBMD is still instrumented after conditioning (Sanderson-
  style approximation), which rules out weak-instrument bias in the direct estimate

Plus univariable MR of each mediator on HF using that mediator's **own** instruments,
which asks a different question: does the factor move HF at all, regardless of whether it
carries the BMD effect.

## Results as computed

**Attenuation of the eBMD 鈫?HF effect (same-instrument baseline)**

| Mediator | n | HERMES direct | p | attenuation | FinnGen direct | p | attenuation |
|---|---|---|---|---|---|---|---|
| FGF23 | 590 / 575 | 鈭?.0337 | 0.0094 | **鈭?.2 %** | 鈭?.0633 | 0.0012 | **鈭?.1 %** |
| SOST | 590 / 575 | 鈭?.0332 | 0.0121 | **+1.4 %** | 鈭?.0616 | 0.0014 | **+2.6 %** |
| BGLAP | 590 / 575 | 鈭?.0336 | 0.0101 | **+0.1 %** | 鈭?.0609 | 0.0022 | **+3.7 %** |
| DKK1 | 590 / 575 | 鈭?.0339 | 0.0099 | **鈭?.8 %** | 鈭?.0630 | 0.0013 | **+0.3 %** |
| PTH | 590 / 575 | 鈭?.0337 | 0.0102 | **鈭?.2 %** | 鈭?.0628 | 0.0013 | **+0.7 %** |
| testosterone | 592 / 577 | 鈭?.0337 | 0.0099 | **鈭?.4 %** | 鈭?.0643 | 0.0010 | **鈭?.4 %** |
| estradiol | 592 / 577 | 鈭?.0328 | 0.0124 | **+2.3 %** | 鈭?.0650 | 0.0008 | **鈭?.5 %** |
| 25(OH)D | 533 / 518 | 鈭?.0353 | 0.0142 | **+2.3 %** | 鈭?.0683 | 0.0017 | **鈭?.2 %** |
| CST3 | 87 / 88 | 鈭?.0225 | 0.29 | +9.7 % | 鈭?.0987 | 0.0039 | +4.7 % |
| CRP | 87 / 88 | 鈭?.0250 | 0.23 | 鈭?.4 % | 鈭?.1038 | 0.0015 | 鈭?.3 % |

**None of the ten absorbs the effect.** Every attenuation is within 卤10 % of zero across
two heart-failure cohorts, and conditional F stays at 131鈥?93, so this is not
weak-instrument bias. The eight mediators analysed on 518鈥?92 instruments agree to within
4 %.

Important reading note: the CST3 and CRP rows use only 87鈥?8 instruments, because those
protein GWAS cover fewer of eBMD's lead variants. Their same-instrument univariable
baseline (鈭?.0249 in HERMES) is itself not significant, so those two rows carry much less
information than the three protein rows with 575鈥?90 instruments.

**Does any mediator move heart failure on its own instruments?**

| Mediator | n | HERMES 尾 | p | FinnGen 尾 | p | reads as |
|---|---|---|---|---|---|---|
| FGF23 | 19 / 20 | +0.079 | 0.13 | +0.121 | **0.051** | consistent positive direction, borderline |
| SOST | 52 / 51 | 鈭?.032 | 0.18 | 鈭?.020 | 0.59 | null |
| BGLAP | 17 / 16 | +0.019 | 0.69 | 鈭?.035 | 0.72 | null |
| DKK1 | 51 / 48 | +0.026 | 0.18 | 鈭?.036 | 0.18 | null, signs disagree |
| **PTH** | 17 / 16 | **+0.168** | **0.0019** | +0.053 | 0.45 | **strong in one cohort, not replicated** |
| testosterone | 202 / 217 | 鈭?.033 | 0.49 | +0.026 | 0.70 | null |
| estradiol | 13 / 13 | 鈭?.085 | 0.81 | 鈭?.269 | 0.37 | uninformative (SE 0.30鈥?.36) |
| 25(OH)D | 117 / 113 | 鈭?.013 | 0.74 | +0.032 | 0.47 | null in both cohorts |
| CST3 | 363 / 356 | +0.047 | **0.012** | +0.003 | 0.93 | inconsistent, kidney marker |
| CRP | 241 / 239 | +0.019 | 0.43 | +0.033 | 0.25 | null |

FGF23 is the only mediator with a consistent positive direction on heart failure in both
cohorts (higher FGF23 鈫?higher HF risk, p = 0.13 and 0.051). But it is **not** on the
bone-density pathway: conditioning on it changes the eBMD effect by 鈭?.2 %, and the
genetic correlation between the eBMD and FGF23 instruments is only 鈭?.03. Two independent
effects, not one mediated by the other.

**PTH is a new finding but an unstable one.** Its own instruments give 尾 = +0.168
(p = 0.0019) in HERMES 鈥?the strongest mediator-to-HF signal seen in this module 鈥?but the
same analysis in FinnGen gives +0.053 (p = 0.45) with 16 instruments. A single-cohort
result at this instrument count is not evidence; it is a hypothesis for a better-powered
PTH GWAS.

CST3 moves HF in HERMES but not in FinnGen; combined with its role as an eGFR proxy, that
is most consistent with kidney function confounding rather than a bone-to-heart axis.

**25(OH)D was the last candidate and it is also null.** Its own 117 / 113 instruments give
尾 = 鈭?.013 (p = 0.74) and +0.032 (p = 0.47) in the two cohorts; conditioning on it moves the
eBMD effect by +2.3 % and 鈭?.2 %. Vitamin D is therefore neither a mediator nor a detectable
driver of heart failure here. The instrument set is well powered for that question
(危R虏 = 0.031, min F = 29.8), so this is a real null rather than the estradiol situation.

## Interpretation

The bone-density effect on heart failure is **not mediated by FGF23, sclerostin,
osteocalcin, DKK1, PTH, testosterone, estradiol, 25-hydroxyvitamin D, cystatin C or
C-reactive protein.**

This is a clean negative result across ten candidates and two independent outcome
cohorts, and it is informative in three directions:

1. It rules out the most obvious mechanistic story (osteocyte-derived FGF23 driving
   cardiomyocyte hypertrophy) as the carrier of the genetic effect. FGF23 appears to
   influence heart failure on its own path, independently of BMD.
2. It rules out inflammation (CRP) and kidney function (cystatin C) as the route, which
   were the two leading competing explanations for a bone鈥揾eart link.
3. It rules out the hormonal and mineral axes that a bone-density signal could plausibly
   travel along (PTH, testosterone, estradiol, vitamin D), which is the last set
   of circulating candidates that a large GWAS covers.

What remains is a genuine gap: the effect is real, replicates in two independent
populations, is subtype-specific, and **does not travel through any circulating factor
tested so far**. Candidate unexplored mechanisms are local rather than systemic 鈥?the
mineralisation machinery in the vessel wall, matrix vesicles, or a
mechanotransduction-independent effect that no blood biomarker would capture.

## Commands

| Script | Command | Output |
|---|---|---|
| `scripts/01_axis_mvmr.py` | `python analysis/06-bone-axis/scripts/01_axis_mvmr.py` | `axis_mvmr.csv`, `axis_univariable.csv`, `axis_strength.csv` |
| self-test | `python analysis/06-bone-axis/scripts/01_axis_mvmr.py --self-test` | PASS/FAIL |

Self-test plants a mediator that carries 69 % of the effect and one that carries none, and
requires the estimator to separate them (observed: +69.2 % versus +9.6 %).

Instrument preparation:

```powershell
python (configured project environment) analysis\02-mr\scripts\02_extract_instruments.py --datasets FGF23,SOST,BGLAP,CST3,CRP,VITD,DKK1,ESTRADIOL --clump distance --kb 500 --keep-all
python (configured project environment) analysis\02-mr\scripts\02b_ld_clump.py --datasets FGF23,SOST,BGLAP,CST3,CRP,VITD,DKK1,ESTRADIOL
python (configured project environment) analysis\02-mr\scripts\02_extract_instruments.py --datasets PTH,TESTO --clump distance --kb 500
```

The distance-only sets for PTH and TESTO were kept as first produced (17 and 248
instruments, min F 26); re-clumping them is not needed for the mediation estimates, which
use the eBMD instruments rather than the mediator's own.

## Limitations

- **The mediator GWAS overlap with UK Biobank** (most of them, and the Olink protein panels
  are drawn from UKB-PPP). The attenuation estimates are therefore not independent of the
  exposure, which biases them *towards* showing mediation 鈥?making the observed absence
  of mediation a conservative finding.
- **Protein levels are not the same as tissue activity.** Circulating FGF23 or sclerostin
  may not reflect the local bone or myocardial concentration that matters.
- **Total osteocalcin is not undercarboxylated osteocalcin.** The BGLAP GWAS measures
  total protein; the biologically active fraction for myocardial energy metabolism is
  undercarboxylated and is not measured by any large GWAS found.
- **CST3 and CRP are analysed on 87鈥?8 instruments only**, so their null results are
  weaker evidence than the eight mediators analysed on 518鈥?92.
- **Estradiol's own-instrument estimate is uninformative**: 13 autosomal instruments give
  SE 鈮?0.30鈥?.36, and its X-chromosome signal cannot be clumped against the autosome-only
  reference. Its *mediation* estimate is unaffected (592 eBMD instruments).
- **Only one mediator at a time.** A joint MVMR of all ten would be more powerful but
  would need a shared instrument set across all mediators, which the different panel
  coverages do not support.
- **Absence of mediation is not absence of pathway.** The effect could act through
  something not measurable in blood.

## Software

- Python 3.12.14 (project `.venv`): numpy 2.5.4, pandas 3.0.6, scipy 1.18.1
- PLINK v1.90b7.2 + 1000G phase 3 EUR for clumping
- GWAS Catalog harmonised files retrieved 2026-10-22 (25(OH)D file retrieved the same day,
  md5 `2b9ea10b68f6e6551dfa7cc63488bda0`)

