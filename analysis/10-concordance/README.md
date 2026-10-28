# 10-concordance 鈥?are the two heart-failure cohorts really telling the same story?

## Purpose

Module `05-hf-followup` reported eBMD 鈫?HF as 尾=鈭?.034 in HERMES (Levin 2022) and
尾=鈭?.064 in FinnGen R12, a twofold difference that was repeatedly described as
"independent replication, stronger in FinnGen". That phrasing assumed the
difference was real. This module tests it.

## What was done

Uses the 588 instruments with complete effects in **both** heart-failure cohorts
(the 600 clumped eBMD instruments minus 12 not present or not usable in one file).

1. **IVW in each cohort** from the same instrument set, so the two estimates are
   directly comparable.
2. **Formal between-cohort test**: z test on the difference of the two IVW
   estimates, using their own standard errors.
3. **Per-instrument concordance**: correlation of the Wald ratio (尾_outcome/尾_exposure)
   across cohorts, Pearson and Spearman.
4. **Precision-weighted regression** of FinnGen on HERMES Wald ratios 鈥?separates a
   proportional difference (slope 鈮?1) from a constant offset (intercept 鈮?0).
5. **Leave-one-instrument-out leverage**: which single variant moves the FinnGen
   estimate most.
6. **Measurement-error correction** (reliability from the ratio SEs). This
   correction **failed** 鈥?see below 鈥?and is reported as a failed diagnostic,
   not as a result.

## Results

| Quantity | Value |
|---|---|
| instruments in both cohorts | 588 |
| HERMES IVW | **鈭?.0345** (se 0.0097) |
| FinnGen IVW | **鈭?.0610** (se 0.0155) |
| ratio | 1.77 |
| **between-cohort difference** | 鈭?.0265 (se 0.0183), z = 鈭?.45, **p = 0.147** |
| slope (FinnGen on HERMES) | 0.609 |
| intercept | 鈭?.0452 |
| Pearson r | 0.425 |
| Spearman 蟻 | 0.352 |
| top-quintile r (117 strongest instruments) | 0.362 |
| largest single-instrument LOCO shift | +0.0053 (rs12251299), 8.7 % relative |

## Interpretation

**The twofold difference is not statistically supported.** The difference is
鈭?.027 with a standard error of 0.018; p = 0.15. Both cohorts agree that the
effect is protective and both are significant on their own, but the estimates are
compatible with a single shared effect. Every statement of the form "the effect is
larger in FinnGen" must be removed or rewritten as "the point estimate is larger
but the difference is not significant".

**The per-instrument agreement is modest** (r = 0.36鈥?.43). This is expected when
most instruments have small exposure effects and therefore noisy Wald ratios in
the outcome cohort; the median instrument contributes little. The top-quintile
subset, where noise should be smaller, does **not** show better agreement (0.362),
which argues that attenuation by noise is not the whole explanation.

**The measurement-error correction failed.** The disattenuated r came out at 1.20,
which is impossible. The reliability model assumes the ratio SE is pure
measurement error in the Wald ratio, but the outcome SE is scaled by 尾_exposure,
which correlates with the ratio itself, so reliability is overestimated and the
correction overshoots. This diagnostic is therefore **discarded**, and no
disattenuated correlation should be quoted.

**The intercept dominates the slope.** slope = 0.61 with intercept = 鈭?.045 means
the cohorts differ mostly by a roughly constant offset rather than by a
proportional scaling of the same per-variant effects. A constant offset is what a
difference in case definition or outcome ascertainment would produce, but with
only two cohorts and p = 0.15 this cannot be separated from chance.

## What this changes in the write-up

- "independently replicated **and stronger**" 鈫?"independently replicated;
  the point estimate is larger but not significantly so (p = 0.15)"
- "effect size nearly doubles (OR 0.938 vs 0.967)" 鈫?keep the numbers, add that
  the difference is within noise
- No claim that per-variant effects are concordant 鈥?they are moderately
  correlated at best (r 鈮?0.4)

## Inputs

- Exposure instruments: `results/tables/02-mr/instruments_eBMD.csv` (600 variants)
- Outcomes: `data/raw/gwas/GCST90162626_HF.h.tsv.gz` (HERMES),
  `data/raw/finngen/finngen_R12_I9_HEARTFAIL.gz` (FinnGen R12)
- Both retrieved 2026-10-21; md5 in `results/tables/download_plan.csv`

## Outputs

- `results/tables/10-concordance/instrument_effect_concordance.csv` 鈥?588 rows,
  per-instrument exposure/outcome effects and Wald ratios in both cohorts
- `results/tables/10-concordance/concordance_summary.csv` 鈥?one row, all scalars above

## Command

```powershell
python (configured project environment) analysis\10-concordance\scripts\01_effect_concordance.py
```

## Dependencies

Python 3.12.14, numpy 2.5.4, pandas 3.0.6, scipy 1.18.1.

## Known limitations

- Only two cohorts, so "offset vs slope" cannot be resolved statistically.
- The instrument set is the clumped eBMD set; instruments absent from either
  outcome file are dropped rather than imputed.
- No correction for winner's curse in the discovery instrument selection.

