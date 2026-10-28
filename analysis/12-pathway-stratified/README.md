# 12-pathway-stratified 鈥?does the diffuse effect sit in any pathway?

## Why this module exists

Module `09-integration` established that the eBMD effect on heart failure is diffuse at
the **locus** level: 589 instruments map to 519 genes, 600 loci contain exactly one shared
causal variant (ABO, PP.H4 = 0.90) and the median PP.H4 is 0.002. Diffuse loci do not
imply diffuse pathways 鈥?a hundred variants in WNT/osteoblast genes could still carry
most of the estimate. This module asks whether any annotated programme carries the effect,
and answers it twice: with pre-specified bone pathways, and with the whole MSigDB
collection.

## What was run

| Script | Question | Output |
|---|---|---|
| `01_pathway_mr.py` | Do canonical bone pathways carry the effect? | `pathway_contributions.csv`, `pathway_group_estimates.csv`, `pathway_heterogeneity.csv`, `pathway_leave_one_out.csv` |
| `02_gene_set_scan.py` | Does **any** MSigDB gene set behave differently from the rest? | `geneset_scan.csv`, `geneset_fwer.csv` |

### 01 鈥?exact decomposition over curated bone pathways

For fixed-effect IVW the estimate is an identity in the per-SNP weights:

```
beta = sum_i w_i bx_i by_i / sum_i w_i bx_i^2,   w_i = 1/sy_i^2
```

so with disjoint instrument groups the effect splits exactly:

```
c_g = sum_{i in g} w_i bx_i by_i / sum_all w_i bx_i^2,    sum_g c_g = beta
```

This is not an approximation, and it is not a re-fit: it states where the estimate comes
from. Seven canonical bone-biology sets (osteoclast resorption, WNT/osteoblast,
matrix/mineralisation, sex hormone, adiposity, height/growth, immune) were applied in a
fixed first-match order; every gene not listed is kept as its own `unassigned` group rather
than dropped. Between-group Cochran Q and a 5,000-draw permutation of the gene labels
(keeping group sizes) test whether the partition explains anything.

### 02 鈥?systematic MSigDB scan

All 9,427 sets from Hallmark, GO Biological Process and Reactome were exported from R
(`msigdbr`) to `data/processed/gene_sets/msigdb_human.csv`. For every set with 鈮?8
instruments, the IVW estimate on that set's instruments is compared with the IVW on the
remaining instruments (disjoint, so SEs add). Multiplicity is handled by Benjamini-Hochberg
across sets **and** by a family-wise permutation that shuffles gene 鈫?instrument assignment
2,000 times and recomputes the largest |z|; the permutation is the honest test because the
MSigDB sets overlap heavily.

## Results

### The canonical bone genes carry essentially none of the effect

HERMES, 592 instruments, fixed-effect 尾 = 鈭?.0336:

| Group | n | share of the effect | own IVW | p |
|---|---|---|---|---|
| **unassigned** | **561** | **+103.6 %** | **鈭?.0389** | **0.0046** |
| wnt_osteoblast | 19 | 鈭?.2 % | +0.0160 | 0.71 |
| sex_hormone | 2 | 鈭?.4 % | 鈥?| 鈥?|
| matrix_mineralisation | 1 | 鈭?.1 % | 鈥?| 鈥?|
| osteoclast_resorption | 4 | 鈭?.4 % | +0.0373 | 0.85 |
| adiposity_bodycomp | 1 | 鈭?.2 % | 鈥?| 鈥?|
| immune_inflammatory | 2 | +2.3 % | 鈥?| 鈥?|
| height_growth | 2 | +2.4 % | 鈥?| 鈥?|

FinnGen gives the same picture (unassigned 546 instruments, +95.5 %, IVW 鈭?.0690,
p = 0.0005; WNT/osteoblast 19 instruments, 鈭?.1 %, IVW +0.0142, p = 0.78).

Two things are worth stating plainly. First, the classic bone loci 鈥?`SOST`, `LRP5`,
`WNT16`, `ESR1`, `RUNX2`, `BGLAP` 鈥?sit on 29 of 592 instruments and their own estimate is
**positive** (opposite sign to the pooled effect), so the effect does not come from them.
Second, between-group Q is 1.63 (df 2, p = 0.44) in HERMES and 3.11 (df 2, p = 0.21) in
FinnGen; the permutation p is 0.45 and 0.26. The curated partition explains nothing.

### The MSigDB scan finds nothing beyond chance

| Cohort | sets 鈮?8 instruments | observed max \|z\| | permutation 95th pct | family-wise p |
|---|---|---|---|---|
| HERMES | 617 (555 GO:BP, 55 Reactome, 7 Hallmark) | 4.52 | 5.70 | **0.40** |
| FinnGen | 624 (561 GO:BP, 56 Reactome, 7 Hallmark) | 3.64 | 4.74 | **0.67** |

Nominal top hits in HERMES are `GOBP_REGULATION_OF_PHOSPHORUS_METABOLISM` (z = +4.52,
BH q = 0.004) and `GOBP_BLOOD_VESSEL_MORPHOGENESIS` (z = +4.15, q = 0.007); in FinnGen they
are `REACTOME_CHROMATIN_ORGANIZATION` (z = 鈭?.64) and `REACTOME_METABOLISM_OF_RNA`
(z = +3.47). **Neither cohort's top set appears in the other's, and neither exceeds the
maximum that random gene labelling produces.** The phosphorus-metabolism hit is the kind of
result the permutation exists to catch: scanning 617 overlapping sets at this instrument
count makes a z of 4.5 an ordinary event.

## Interpretation

The effect is diffuse in pathways as well as in loci. Three statements are now supported by
measurement rather than by inference from colocalisation:

1. The canonical bone pathway genes carry ~5 % of the estimate, and their own estimate has
   the opposite sign. The instrument set is not "a bone-axis instrument set" in the
   textbook sense.
2. No Hallmark, GO:BP or Reactome set behaves differently from the rest of the genome.
3. This is consistent with the failure of all ten circulating mediators (module 06) and
   explains it: there is no single pathway to absorb, because no annotated pathway is
   over-represented among the effect-carrying variants.

The practical consequence is a negative prediction worth writing down: any further
"eBMD 鈫?X pathway 鈫?heart failure" hypothesis built on nearest-gene annotation starts from
a partition whose largest annotated component accounts for 4 % of the effect. The data
cannot support that kind of claim at this scale.

## Limitations

- **`nearest_genes` is proximity annotation, not causal gene assignment.** 561 of 592
  instruments fall outside the curated lists and most of their nearest genes are
  unannotated or unrelated open-reading frames; the gene-set scan is limited by that.
- **Only 7 Hallmark sets clear the 8-instrument floor.** Hallmark genes are largely
  housekeeping and immune programmes, which is why the GO:BP collection (555 sets) carries
  the actual scan. A negative Hallmark result alone would have been uninformative.
- **The contribution decomposition is fixed-effect.** Random-effects IVW is the headline
  estimate elsewhere in the project; the decomposition uses fixed-effect weights so that
  the shares sum exactly. The two differ by less than 0.1 % here (尾_FE 鈭?.0336 vs random
  鈭?.0336).
- **Gene sets overlap**, so BH q-values across 600+ sets are optimistic; that is why the
  family-wise permutation, not the q-value, is reported as the test.
- **Power.** With 8鈥?0 instruments per set, only a large difference (z 鈮?5 under this null)
  is detectable. A pathway carrying 10 % of a 尾 of 鈭?.064 is not.

## Commands

```powershell
# export the gene sets once (R)
Rscript -e "library(msigdbr); ..."   # see analysis/12-pathway-stratified/README.md above

python (configured project environment) analysis\12-pathway-stratified\scripts\01_pathway_mr.py --n-perm 5000
python (configured project environment) analysis\12-pathway-stratified\scripts\02_gene_set_scan.py --n-perm 2000

# self-tests
python (configured project environment) analysis\12-pathway-stratified\scripts\01_pathway_mr.py --self-test
python (configured project environment) analysis\12-pathway-stratified\scripts\02_gene_set_scan.py --self-test
```

Self-tests: `01` plants two groups (one carrying 98.5 % of the effect) and requires the
decomposition to sum to 尾 and Q to separate them (observed Q = 9.48, p = 0.002); `02`
plants a 60-instrument set with a doubled effect and requires the scan to rank it first
(observed z = 鈭?.15, p = 0.032).

## Software

- Python 3.12.14 (project `.venv`): numpy 2.5.4, pandas 3.0.6, scipy 1.18.1
- R 4.6.1 with `msigdbr` for the gene-set export (Hallmark, GO:BP, Reactome)

