# 13-robustness 鈥?three tests that decide how the finding may be worded

> **Panel A below has been superseded. Read `analysis/15-modality/README.md` first.**
> This module tested one DXA file and found a significant estimate of the opposite sign.
> Module 15 added the two other DXA femoral-neck files that exist: they give +0.039 (n=47,205)
> and +0.025 (n=32,017), neither significant, and in HERMES the third is **negative**
> (鈭?.039). The DXA pool is +0.047 (p=0.019) in FinnGen but +0.002 (p=0.92, I虏=0.60) in HERMES.
> The correct statement is therefore **"DXA does not show the effect and the one positive
> estimate did not replicate"**, not "DXA shows the opposite". The wording "heel-ultrasound
> bone density" still stands 鈥?but on the absence of positive DXA evidence, not on evidence of
> an opposite DXA effect. Panels B and C of this module are unaffected.

## Why this module exists

Every genetic estimate in modules 02鈥?2 rests on two unstated assumptions:

1. **the exposure is bone density**, when in fact it is *heel ultrasound* eBMD from UK
   Biobank (GCST90029004, n=445,921) 鈥?a speed-of-sound phenotype, not DXA mineral content;
2. **the outcome is reached through the myocardium**, which was inferred by elimination
   (not coronary, not circulating, not body size, not any annotated pathway) rather than
   measured.

Both are testable with public data. A third question 鈥?whether the instrument set has
generic pleiotropy, which would undermine "endpoint-specific" 鈥?is testable with negative
controls. All three decision rules were written into the script headers before the runs.

| Script | Question | Decision rule |
|---|---|---|
| `01_dxa_exposure.py` | Does DXA femoral neck BMD give the same effect? | direction agrees + CI excludes 0 鈫?the effect belongs to bone density; null with CI excluding 鈭?.064 鈫?the finding is specific to the heel-ultrasound phenotype |
| `02_lv_mass.py` | Does eBMD move left ventricular mass? | significant negative 鈫?cardiac remodelling is measured evidence; null with a tight CI 鈫?the axis is not concentric hypertrophy |
| `03_negative_controls.py` | Is the instrument set generically pleiotropic? | any Bonferroni-significant negative control 鈫?withdraw "endpoint-specific" |

## A. The two modalities disagree 鈥?this is the module's headline

**Input.** DXA femoral neck BMD, GCST90446641 (n=79,940 UK Biobank; PubMed 39349818).
6,548 variants at p<5e-8 鈫?PLINK clumped (r虏<0.001, 10 Mb, 1000G EUR) to **88 instruments**,
危R虏 = 0.0737, min F = 29.7. Canonical BMD genes are present in the set 鈥?`WNT16`, `WLS`,
`SOST`, `HDAC5`, `SOX6`, `MEF2C`, `COLEC10`, `CCDC170` 鈥?so the file is genuine femoral neck
BMD and not a mislabelled trait.

| Outcome | 尾 (95 % CI) | p | power for \|尾\|=0.064 |
|---|---|---|---|
| **DXA femoral neck BMD 鈫?FinnGen heart failure** | **+0.0678 (+0.0052, +0.1303)** | **0.034** | 52 % |
| DXA femoral neck BMD 鈫?HERMES heart failure | +0.0321 (鈭?.0130, +0.0772) | 0.163 | 79 % |
| 鈥?heel-ultrasound eBMD 鈫?FinnGen (for contrast) | 鈭?.0640 (鈭?.1005, 鈭?.0276) | 0.00057 | 93 % |

**The two modalities point in opposite directions, and the FinnGen CIs do not overlap.**

The obvious explanations were tested and none of them removes the disagreement:

| Check | Result | Reading |
|---|---|---|
| Allele coding | FNBMD instruments 鈫?eBMD **+0.539** (se 0.062), Egger intercept p=0.51; eBMD instruments 鈫?FNBMD **+0.469** (se 0.021), p=0.068 | both files describe real, positively correlated BMD; the sign is not a coding error |
| Single-locus dependence | leave-one-out range **+0.056 to +0.078**, never crossing zero; dropping the dominant locus (rs11551167, 50 % of the weight) gives +0.0686 (p=0.043) | not one variant |
| Body size | MVMR: +BMI +0.064 (p=0.048); +lean mass +0.052 (p=0.089); +height +0.073 (p=0.028); +BMI+height +0.067 (p=0.074); BMI-near loci removed (12) +0.065 (p=0.065) | not adiposity, not stature |
| Generic pleiotropy | 0 of 6 negative controls pass Bonferroni; 1 nominal (appendicitis +0.059, p=0.024; 0.3 expected by chance) | no gross instrument failure |
| Subtypes | HF+CHD +0.033 (p=0.49); HF+obesity +0.049 (p=0.22) | no subtype specificity either way |

`dxa_positional_concordance.csv` reports a low effect correlation at positionally matched
loci (r = 0.078 over 28 allele-alignable pairs). That check is weak by construction 鈥?the
600 eBMD instruments are spaced ~5 Mb apart, so a "nearest neighbour within 500 kb" is
usually a different variant 鈥?and it should not be read as evidence against the two traits
being correlated. The reciprocal-prediction estimates above are the usable concordance
measure, and they are strongly positive.

**What this means.** The protective effect is established for the **heel-ultrasound**
phenotype and does not transfer to DXA areal BMD at the femoral neck. Two readings remain
open, and the data in hand cannot separate them:

- **biological**: quantitative ultrasound reflects trabecular architecture and elasticity,
  which areal BMD does not capture; the causal property would then be bone *quality*, not
  mineral content;
- **technical**: the DXA summary statistics come from a 2024 methods paper (machine-learning
  assisted GWAS inference), and no second large DXA femoral neck GWAS with summary
  statistics exists to check them against.

Either way, the correct wording is now "heel-ultrasound bone density", and the DXA
disagreement belongs in the abstract, not in a supplementary caveat.

## B. Left ventricular mass: small, negative, directionally consistent

| Trait | n instruments | raw 尾 | SD-standardised 尾 (95 % CI) | p | power for 0.04 SD |
|---|---|---|---|---|---|
| Left ventricular mass (g) | 591 | 鈭?.969 g | **鈭?.041 SD (鈭?.074, 鈭?.008)** | **0.015** | 66 % |
| LV mass indexed to body surface area | 591 | 鈭?.099 | 鈭?.016 SD (鈭?.049, +0.017) | 0.34 | 66 % |

Higher genetically predicted eBMD goes with **lower** left ventricular mass, at about
0.04 SD per SD 鈥?a small structural effect, nominally significant in the unindexed trait
and attenuated (still negative) once body surface area is divided out. Both files are UK
Biobank CMR, so they overlap the exposure; overlap inflates false positives, which makes
the correct reading "consistent with, not proof of" a myocardial route.

**The LV-mass 鈫?heart-failure direction is not estimable**: the LV mass GWAS yields only 9
independent instruments, below the module's floor of 10.

## C. Negative controls: the endpoint specificity holds

Six FinnGen R12 endpoints pre-specified as unrelated to bone density, all well powered
(26,894鈥?3,886 cases):

| Endpoint | heel eBMD 尾 (p) | DXA neck BMD 尾 (p) |
|---|---|---|
| acute appendicitis | +0.034 (0.052) | +0.059 (**0.024**) |
| senile cataract | 鈭?.005 (0.77) | +0.015 (0.53) |
| atopic dermatitis | 鈭?.005 (0.81) | 鈭?.036 (0.21) |
| migraine | +0.013 (0.53) | +0.005 (0.87) |
| tonsillectomy / adenotomy | +0.009 (0.64) | +0.008 (0.77) |
| acute sinusitis | 鈭?.025 (0.20) | 鈭?.012 (0.75) |

**heel eBMD: 0 of 6 nominal (minimum p = 0.052), 0 Bonferroni, median power 92 %.** The
instrument set is not generically pleiotropic, so "endpoint-specific" survives.

**DXA: 1 of 6 nominal** (appendicitis, p = 0.024; binomial upper p = 0.26 鈥?ordinary for six
tests), 0 Bonferroni, median power 64 %. The DXA set is not obviously broken either, which
is precisely why the disagreement in panel A cannot be dismissed as an invalid instrument
set.

Worth recording: eBMD's appendicitis estimate (+0.034, p=0.052) sits exactly at the
conventional threshold and points the same way as its heart-failure estimate. It does not
survive correction, but it is the one number in this table that could become interesting
with more power.

## What this module establishes

1. **The exposure must be named.** The finding is "heel-ultrasound BMD and heart failure",
   not "bone density and heart failure". A DXA femoral neck GWAS of adequate size gives a
   significantly opposite estimate, and every diagnostic available here leaves that
   disagreement standing.
2. **"Endpoint-specific" is safe.** Six unrelated endpoints are null for the eBMD set at
   92 % median power.
3. **The myocardial claim improves from inference to measurement, weakly.** LV mass falls
   by 0.04 SD per SD of eBMD; it is small, UKB-overlapped, and not the concentric
   hypertrophy pattern.
4. **Not testable with public data**: the clean "heart failure without coronary disease"
   contrast (FinnGen R12 has only HF, HF+CHD, HF+hypertrophic cardiomyopathy [561 cases]
   and HF+overweight 鈥?there is no HF-excluding-CHD endpoint), and LV mass 鈫?heart failure
   (9 instruments).

## Outputs

| File | Contents |
|---|---|
| `dxa_exposure.csv` | DXA BMD 鈫?HF and 鈫?FinnGen HF with power |
| `dxa_concordance.csv` | rsID-shared instrument comparison |
| `dxa_positional_concordance.csv`, `dxa_matched_loci.csv` | 卤500 kb nearest-neighbour comparison |
| `dxa_instrument_genes.csv` | nearest gene for each DXA instrument |
| `dxa_diagnostics_mvmr.csv` | DXA BMD adjusted for BMI / lean mass / height, and with BMI-near loci dropped |
| `lv_mass.csv` | eBMD 鈫?LV mass and LV mass/BSA, raw and SD-standardised |
| `negative_controls.csv`, `negative_summary.csv` | six endpoints 脳 two exposures, Bonferroni and binomial summaries |

Figure: `figures/fig8_robustness.png`.

## Commands

```powershell
python (configured project environment) analysis\scripts\download_plan.py --only FNBMD,LVM,LVMI,FN_NEG_APPENDACUT,FN_NEG_CATARACT,FN_NEG_ATOPIC,FN_NEG_MIGRAINE,FN_NEG_TONSILL,FN_NEG_SINUSITIS --chunk-mb 4 --connections 8
python (configured project environment) analysis\02-mr\scripts\02_extract_instruments.py --datasets FNBMD,LVM --clump distance --kb 500 --keep-all
python (configured project environment) analysis\02-mr\scripts\02b_ld_clump.py --datasets FNBMD,LVM

python (configured project environment) analysis\13-robustness\scripts\01_dxa_exposure.py
python (configured project environment) analysis\13-robustness\scripts\02_lv_mass.py
python (configured project environment) analysis\13-robustness\scripts\03_negative_controls.py
python (configured project environment) analysis\13-robustness\scripts\04_dxa_diagnostics.py
python (configured project environment) analysis\scripts\make_figure8.py
```

Self-tests (`--self-test` on each script) check: allele alignment flips signs correctly;
a planted 鈭?.064 is recovered and a null is not; a planted effect inside a 592-SNP set is
detected while six null endpoints stay null.

## Limitations

- **The DXA file is a secondary deposit from a methods paper.** Its instruments behave
  correctly in every check run here, but no independent DXA femoral neck GWAS exists at
  comparable size to confirm the reversal. This is the single biggest open question in the
  project.
- **Both LV mass files are UK Biobank**, so they share participants with the exposure.
  Overlap biases towards false positives; the negative LV-mass-indexed estimate should be
  weighted accordingly.
- **Negative controls cannot detect pleiotropy that acts through correlated traits.** They
  exclude gross instrument failure, not subtle confounding by, for example, general health.
- **Power for the DXA estimates is 52鈥?9 %** for the effect size of interest, so "the DXA
  effect is positive" is itself only moderately supported; what is solid is that it is
  *not* the same as the eBMD effect.

## Software

- Python 3.12.14 (project `.venv`): numpy 2.5.4, pandas 3.0.6, scipy 1.18.1
- PLINK v1.90b7.2 + 1000G EUR for clumping
- FinnGen R12 endpoint manifest read from `data/processed/finngen_R12_manifest.tsv`
  (2,469 endpoints; case/control counts used for power)

