# 09-integration 鈥?joining the genetic, cell-type and colocalisation layers

## Purpose

Modules 02鈥?8 produced three things that never touched each other:

- 600 LD-clumped eBMD instruments, with a heart-failure effect of OR 0.938 per SD in FinnGen
- 9 circulating candidate mediators, every one of which failed to absorb the effect
- a single-cell map of calcified carotid plaque in which the osteogenic transcription-factor
  programme sits in macrophages while bone-pathway genes sit in fibroblasts

Nothing connected them. A variant is not a gene; a gene's compartment is not a mechanism;
and a shared locus is not a shared causal variant. This module makes those three joins
explicit and tests each one.

## Step 1 鈥?instruments to genes (`01_instrument_genes.py`)

FinnGen summary statistics carry a `nearest_genes` column, so assignment needs no external
annotation and no network: scan the FinnGen HF file once and read the gene off the rows
whose rsID is an eBMD instrument.

**Result: 589 of 600 instruments (98.2 %) matched**, mapping to 519 distinct nearest genes.
Known bone loci come back out of the annotation unprompted 鈥?CPED1, ESR1, SOST, SOX6,
TNFSF11 (RANKL) 鈥?which is the check that the join is doing what it claims.

## Step 2 鈥?genes to compartments (`02_gene_to_celltype.py`)

428 of the 589 genes (72.7 %) are measurable in the GSE159677 matrix. For each, mean
expression per cell type; each gene's *preferred compartment* is the cell type where it is
highest, counted only if detected in 鈮?0 % of some compartment (194 genes).

The null is matched, not random: instrument genes sit at disease loci, which are gene-dense
and often long, so comparing against all genes would confound "instrument gene" with
"expression level". The background is 1,000 random gene sets drawn within expression
deciles.

| Compartment | observed share | matched null (mean) | fold | p |
|---|---|---|---|---|
| **fibroblast** | 0.263 | 0.200 | **1.31** | **0.016** |
| **smooth muscle** | 0.217 | 0.166 | **1.31** | **0.039** |
| endothelial | 0.217 | 0.234 | 0.93 | 0.70 |
| **macrophage** | 0.217 | 0.247 | **0.88** | 0.84 |
| T cell | 0.052 | 0.077 | 0.67 | 0.91 |
| B cell / plasma | 0.036 | 0.070 | 0.52 | 0.97 |
| mast | 0.000 | 0.007 | 0.00 | 1.00 |

**The genes the BMD instruments point at are enriched in fibroblasts and smooth muscle 鈥?the two mesenchymal compartments 鈥?and are, if anything, depleted from macrophages**,
which is where the osteogenic transcription-factor programme sits. Those are two different
cell populations carrying two different aspects of "bone-like" behaviour, and the genetic
signal follows the mesenchymal one.

## Step 3 鈥?colocalisation (`03_coloc.py`)

A shared nearest gene is not a shared mechanism. coloc-abf (Giambartolomei 2014) with at
most one causal variant per trait per locus was run on 卤500 kb around each of the 600
instruments: 1,871,584 variants shared between eBMD and FinnGen HF after allele alignment
and palindromic filtering. Wakefield ABF with sd_prior 0.15 (quantitative) and 0.2 (case-
control), sd_Y = 1/鈭?s(1鈭抯)) with s = 37,653/500,348; coloc's default priors p1=p2=1e-4,
p12=1e-5.

| Class | Loci |
|---|---|
| only one trait associated (PP.H1 or H2 > 0.5) | 583 / 600 |
| both associated, **different** causal variants (PP.H3 鈮?0.8) | 6 |
| **shared** causal variant (PP.H4 鈮?0.8) | **1** |
| median PP.H4 across all loci | **0.002** |

The six PP.H3 loci are SMG6, FHL3, GHR, ZFHX3, FAF1 and NFIA 鈥?both traits associate, but
through distinct variants, so they are not shared mechanisms.

**The single colocalising locus is ABO** (rs635634, chr9:133,279 Mb, PP.H4 = 0.982). Both
traits are strongly associated there in the same direction: eBMD 尾 = +0.0283 (z = 11.6),
FinnGen HF 尾 = 鈭?.0484 (z = 鈭?.8).

## Step 4 鈥?does ABO carry the effect? (`04_abo_exclusion.py`)

ABO is the most pleiotropic locus in the genome 鈥?blood-group glycosylation shifts von
Willebrand factor, factor VIII, E-selectin, ICAM-1 and LDL. A shared causal variant there
is at least as likely to be horizontal pleiotropy as a bone-to-heart mechanism, so the
decisive question is whether removing it removes the effect.

| Cohort | Subset | n | 尾 (95 % CI) | p |
|---|---|---|---|---|
| HERMES | all instruments | 592 | 鈭?.0336 (鈭?.059, 鈭?.008) | 0.010 |
| HERMES | ABO region removed | 591 | 鈭?.0329 (鈭?.058, 鈭?.008) | 0.011 |
| HERMES | PP.H3 loci removed | 586 | 鈭?.0393 (鈭?.065, 鈭?.014) | 0.0023 |
| HERMES | both removed | 585 | 鈭?.0386 (鈭?.064, 鈭?.014) | 0.0025 |
| FinnGen | all instruments | 577 | 鈭?.0640 (鈭?.101, 鈭?.028) | 0.0006 |
| FinnGen | ABO region removed | 576 | 鈭?.0639 (鈭?.100, 鈭?.027) | 0.0006 |
| FinnGen | PP.H3 loci removed | 571 | 鈭?.0664 (鈭?.103, 鈭?.030) | 0.0003 |
| FinnGen | both removed | 570 | 鈭?.0663 (鈭?.103, 鈭?.030) | 0.0003 |

**ABO carries essentially none of the effect**: only one instrument sits within 1 Mb of it,
and removing it moves the estimate by 0.0001 in FinnGen and 0.0007 in HERMES. Removing the
PP.H3 loci makes the estimate *stronger* in both cohorts.

## What this module establishes

The bone-density effect on heart failure is **diffuse and polygenic, not localisable**:

- 519 genes across the genome, no locus dominant
- median PP.H4 = 0.002: at almost no instrument locus is there a shared causal variant
- the one locus that does colocalise (ABO) contributes 鈮?0 to the effect
- leave-one-locus-out over the strongest variants moves the estimate by 鈮?0.011
- the genes involved are expressed in the mesenchymal compartments, not the macrophage
  compartment that carries the osteogenic transcription programme

This is the mechanistic explanation for why ten circulating mediators all failed: there is
no single pathway to absorb, because the effect is the sum of hundreds of small ones.

## Commands

| Script | Command | Output |
|---|---|---|
| `scripts/01_instrument_genes.py` | `python 鈥?01_instrument_genes.py` | `instrument_genes.csv` |
| `scripts/02_gene_to_celltype.py` | `python 鈥?02_gene_to_celltype.py --n-null 1000` | `gene_celltype_expression.csv`, `gene_preferred_celltype.csv`, `celltype_enrichment.csv`, `expression_summary.csv` |
| `scripts/03_coloc.py` | `python 鈥?03_coloc.py` | `coloc_loci.csv`, `coloc_summary.csv` |
| `scripts/04_abo_exclusion.py` | `python 鈥?04_abo_exclusion.py` | `abo_exclusion.csv` |
| self-tests | `--self-test` on each | PASS/FAIL |

Self-tests: the gene scanner returns only requested rsIDs; the compartment scorer assigns
100 % of genes planted in a compartment to that compartment; coloc gives PP.H4 = 1.000 for
a shared variant and PP.H3 = 1.000 for two distinct ones; the exclusion routine recovers
the clean estimate after removing planted spurious instruments.

## Limitations

- **`nearest_genes` is a proximity annotation, not a causal gene.** A lead variant may
  regulate a gene megabases away; the assignment is a hypothesis generator, not a result.
- **27 % of instrument genes are not measurable** in the plaque matrix, so the compartment
  enrichment describes the measurable subset and could be biased toward well-expressed,
  longer genes.
- **Single causal variant per locus.** coloc-abf assumes one causal variant per trait per
  window; with 卤500 kb windows containing multiple independent signals this is violated,
  and PP.H4 can be diluted (making the "1 of 600" an underestimate rather than an
  overestimate).
- **ABO cannot be resolved here.** Whether the ABO signal is blood-group pleiotropy or a
  genuine shared mechanism needs stratified analysis by blood group, which summary
  statistics cannot provide.
- **Three donors** in the single-cell data, so the compartment enrichment is a claim about
  these three plaques.

## Software

- Python 3.12.14 (project `.venv`): numpy 2.5.4, pandas 3.0.6, scipy 1.18.1,
  scikit-learn 1.9.1
- coloc-abf reimplemented from Giambartolomei et al. 2014 (PLoS Genet 10:e1004383)

