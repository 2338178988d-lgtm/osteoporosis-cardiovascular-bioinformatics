"""Build the evidence-backed JBMR supplementary material package.

This script only copies/derives from project result tables and recorded metadata;
it does not rerun the scientific analyses.
"""
from __future__ import annotations

import csv
import hashlib
import json
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    KeepTogether,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "supplementary"
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "results" / "reports"
ENV = REPORTS / "鐜璁板綍_environment.txt"


def read_csv(rel: str) -> list[dict[str, str]]:
    p = ROOT / rel
    with p.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def copy_table(src: str, dest: str):
    p = ROOT / src
    q = OUT / "tables" / dest
    q.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(p, q)


def tidy(value: object) -> str:
    return str(value).replace("\n", " ").replace("\r", " ")


# The current evidence-backed mapping. A row is retained even when an accession is
# a study record rather than a separate paper; that distinction is explicit.
REFERENCE_ROWS = [
    ["GCST90029004", "heel-ultrasound bone density", "445,921 in the harmonised analysis file; 583,314 in GWAS Catalog/OpenGWAS discovery metadata", "29892013", "Loh PR et al. Mixed-model association for biobank-scale datasets. Nat Genet. 2018;50:906-908.", "10.1038/s41588-018-0144-6", "analysis file n column; GWAS Catalog study page; OpenGWAS metadata"],
    ["GCST90025982", "heel-ultrasound BMD T-score", "445,855", "34226706", "Barton AR et al. Whole-exome imputation within UK Biobank powers rare coding variant association and fine-mapping analyses. Nat Genet. 2021;53:1260-1269.", "10.1038/s41588-021-00892-1", "GWAS Catalog study page"],
    ["GCST006979", "heel-ultrasound BMD", "426,824", "30598549", "Morris JA et al. An atlas of genetic influences on osteoporosis in humans and mice. Nat Genet. 2019;51:258-266.", "10.1038/s41588-018-0302-x", "GWAS Catalog study page"],
    ["GCST90446641; GCST90446628", "DXA femoral-neck BMD", "79,940; 47,205", "39349818", "Miao J et al. Valid inference for machine learning-assisted genome-wide association studies. Nat Genet. 2024;56:2361-2369.", "10.1038/s41588-024-01934-0", "GWAS Catalog study pages"],
    ["GCST90568458", "DXA femoral-neck BMD", "32,017", "41206123", "Qian Y et al. Genome-wide Association Studies of over 30,000 Samples with Bone Mineral Density at Multiple Skeletal Sites and Its Clinical Relevance. Genomics Proteomics Bioinformatics. 2025;23(5).", "10.1093/gpbjnl/qzaf097", "GWAS Catalog study page"],
    ["GCST90162626", "heart failure (HERMES)", "1,366,492 discovery; 328,071 replication", "36376295", "Levin MG et al. Genome-wide association and multi-trait analyses characterize the common genetic architecture of heart failure. Nat Commun. 2022;13:6914.", "10.1038/s41467-022-34216-6", "GWAS Catalog study page; PubMed"],
    ["GCST009541", "heart failure replication", "47,309 cases; 930,014 controls", "31919418", "Shah S et al. Genome-wide association and Mendelian randomisation analysis provide insights into the pathogenesis of heart failure. Nat Commun. 2020;11:163.", "10.1038/s41467-019-13690-5", "GWAS Catalog study page; PubMed"],
    ["GCST90018806; GCST90018877; GCST90018887; GCST90018586; GCST90018667", "Sakaue cross-population HF/MI/OP/East-Asian traits", "accession-specific counts in Table S1", "34594039", "Sakaue S et al. A cross-population atlas of genetic associations for 220 human phenotypes. Nat Genet. 2021;53:1415-1424.", "10.1038/s41588-021-00931-x", "GWAS Catalog study pages; PubMed"],
    ["GCST90132315", "coronary artery disease", "1,378,170 total discovery label", "36474045", "Aragam KG et al. Discovery and systematic characterization of risk variants and genes for coronary artery disease in over a million participants. Nat Genet. 2022;54:1803-1815.", "10.1038/s41588-022-01233-6", "GWAS Catalog study page; PubMed"],
    ["GCST90018877", "myocardial infarction", "461,823 European + 161,206 East Asian", "34594039", "Sakaue S et al. A cross-population atlas of genetic associations for 220 human phenotypes. Nat Genet. 2021;53:1415-1424.", "10.1038/s41588-021-00931-x", "GWAS Catalog study page"],
    ["GCST90104534", "stroke", "1,308,460 European discovery label", "36180795", "Mishra A et al. Stroke genetics informs drug discovery and risk prediction across ancestries. Nature. 2022;611:115-123.", "10.1038/s41586-022-05165-3", "GWAS Catalog study page; PubMed"],
    ["GCST90244711", "left ventricular mass", "43,230 total discovery label", "36944631", "Khurshid S et al. Clinical and genetic associations of deep learning-derived cardiac magnetic resonance-based left ventricular mass. Nat Commun. 2023;14:1558.", "10.1038/s41467-023-37173-w", "GWAS Catalog study page; PubMed"],
    ["GCST90435260", "left ventricular mass indexed to body surface area", "36,083", "39966646", "Tadros R et al. Large-scale genome-wide association analyses identify novel genetic loci and mechanisms in hypertrophic cardiomyopathy. Nat Genet. 2025;57:530-538.", "10.1038/s41588-025-02087-4", "GWAS Catalog study page; PubMed"],
    ["GCST90000025", "appendicular lean mass", "450,243", "33097823", "Pei YF et al. The genetic architecture of appendicular lean mass characterized by association analysis in the UK Biobank study. Commun Biol. 2020;3:608.", "10.1038/s42003-020-01334-0", "GWAS Catalog study page; PubMed"],
    ["GCST90029007; GCST90029008", "BMI and height", "532,396; 673,878", "29892013", "Loh PR et al. Mixed-model association for biobank-scale datasets. Nat Genet. 2018;50:906-908.", "10.1038/s41588-018-0144-6", "GWAS Catalog study pages"],
    ["GCST90000618", "serum 25-hydroxyvitamin D", "496,946", "32242144", "Revez JA et al. Genome-wide association study identifies 143 loci associated with 25 hydroxyvitamin D concentration. Nat Commun. 2020;11:1647.", "10.1038/s41467-020-15421-7", "GWAS Catalog study page; PubMed"],
    ["GCST90012114; GCST90012105", "testosterone and estradiol", "425,097; 206,927", "32042192", "Ruth KS et al. Using human genetics to understand the disease impacts of testosterone in men and women. Nat Med. 2020;26:252-258.", "10.1038/s41591-020-0751-5", "GWAS Catalog study pages; PubMed"],
    ["GCST90025945; GCST90025959", "cystatin C and C-reactive protein", "437,846; 436,939", "34226706", "Barton AR et al. Whole-exome imputation within UK Biobank powers rare coding variant association and fine-mapping analyses. Nat Genet. 2021;53:1260-1269.", "10.1038/s41588-021-00892-1", "GWAS Catalog study pages"],
    ["GCST90469223; GCST90470711; GCST90468440; GCST90468997; GCST90470372", "FGF23, SOST, BGLAP, DKK1 and PTH", "47,745 each", "39789286", "Loya H et al. A scalable variational inference approach for increased mixed-model association power. Nat Genet. 2025;57:461-468.", "10.1038/s41588-024-02044-7", "GWAS Catalog study pages; PubMed"],
    ["GCST90668010; GCST90668011", "East-Asian HFrEF/HFpEF", "201,831", "41184235", "Enzan N et al. Genome-wide analysis of heart failure yields insights into disease heterogeneity and enables prognostic prediction in the Japanese population. Nat Commun. 2025;16:9680.", "10.1038/s41467-025-64659-6", "GWAS Catalog study pages; PubMed"],
    ["GCST90726619; GCST90726620", "European HFpEF/HFrEF", "756,189; 940,805", "40038546", "Henry A et al. Genome-wide association study meta-analysis provides insights into the etiology of heart failure and its subtypes. Nat Genet. 2025;57:815-828.", "10.1038/s41588-024-02064-3", "GWAS Catalog study pages; PubMed"],
    ["GCST90310293", "calcific aortic valve stenosis", "956,682 total discovery label", "38494474", "Th茅riault S et al. Integrative genomic analyses identify candidate causal genes for calcific aortic valve stenosis involving tissue-specific regulation. Nat Commun. 2024;15:2407.", "10.1038/s41467-024-46639-4", "GWAS Catalog study page; PubMed"],
    ["GCST90475970; GCST90475949; GCST90475951", "VA Million Veteran Program AF/cardiomyopathy", "accession-specific counts in Table S1", "39024449", "Verma A et al. Diversity and scale: Genetic architecture of 2068 traits in the VA Million Veteran Program. Science. 2024;385:eadj1182.", "10.1126/science.adj1182", "GWAS Catalog study pages; PubMed"],
    ["GCST90624411", "cross-population atrial fibrillation", "252,438", "40645996", "Yuan S et al. Cross-population GWAS and proteomics improve risk prediction and reveal mechanisms in atrial fibrillation. Nat Commun. 2025;16:6426.", "10.1038/s41467-025-61720-2", "GWAS Catalog study page; PubMed"],
    ["GCST90455657", "large European heart failure", "168,363 cases", "40195560", "Lee DSM et al. Common-variant and rare-variant genetic architecture of heart failure across the allele-frequency spectrum. Nat Genet. 2025;57:829-838.", "10.1038/s41588-025-02140-2", "GWAS Catalog study page; PubMed"],
    ["GCST90473571; GCST90473575; GCST90473562", "non-Finnish European AF/HF/CM", "accession-specific counts in Table S1", "40770095", "UK Biobank Whole-Genome Sequencing Consortium. Whole-genome sequencing of 490,640 UK Biobank participants. Nature. 2025;645:692-701.", "10.1038/s41586-025-09272-9", "GWAS Catalog study pages; PubMed"],
    ["GSE100927", "peripheral arterial bulk transcriptome", "104", "29500419", "Steenman M et al. Identification of genomic differences among peripheral arterial beds in atherosclerotic and healthy arteries. Sci Rep. 2018;8:3940.", "10.1038/s41598-018-22292-y", "NCBI GEO series page; PubMed"],
    ["GSE43292", "carotid atheroma bulk transcriptome", "64", "23660665", "Ayari H, Bricca G. Identification of two genes potentially associated in iron-heme homeostasis in human carotid plaque using microarray analysis. J Biosci. 2013;38:311-315.", "10.1007/s12038-013-9310-2", "NCBI GEO series page; PubMed"],
    ["GSE28829", "early/advanced carotid plaque", "29", "22388324", "D枚ring Y et al. Auto-antigenic protein-DNA complexes stimulate plasmacytoid dendritic cells to promote atherosclerosis. Circulation. 2012;125:1673-1683.", "10.1161/CIRCULATIONAHA.111.046755", "NCBI GEO series page; PubMed"],
    ["GSE35956", "osteoporosis-derived hMSC", "10", "23028809", "Benisch P et al. The transcriptional profile of mesenchymal stem cell populations in primary osteoporosis is distinct and shows overexpression of osteogenic inhibitors. PLoS One. 2012;7:e45142.", "10.1371/journal.pone.0045142", "NCBI GEO series page; PubMed"],
    ["GSE56815", "bone-side blood monocytes, HG-U133A", "80", "29330445", "Liu Y et al. A novel approach for correction of crosstalk effects in pathway analysis and its application in osteoporosis research. Sci Rep. 2018;8:668.", "10.1038/s41598-018-19196-2", "NCBI GEO series page; PubMed"],
    ["GSE56814", "bone-side blood monocytes, HuEx-1_0-st", "73", "30056508", "Liu Y et al. Transcriptomic Data Identified Key Transcription Factors for Osteoporosis in Caucasian Women. Calcif Tissue Int. 2018;103:581-588.", "10.1007/s00223-018-0457-6", "NCBI GEO series page; PubMed"],
    ["GSE7158", "peak-bone-mass blood monocytes", "26", "19223260", "Lei S et al. An in vivo genome wide gene expression study of circulating monocytes suggested GBP1, STAT1 and CXCL10 as novel risk genes for the differentiation of peak bone mass. Bone. 2009;44:1010-1014.", "10.1016/j.bone.2008.05.016", "NCBI GEO series page; PubMed"],
    ["GSE7429", "osteoporosis blood B cells", "20", "18433299", "Ma M et al. In vivo genome-wide expression study on human circulating B cells suggests a novel ESR1 and MAPK3 network for postmenopausal osteoporosis. J Bone Miner Res. 2008;23:644-654.", "10.1359/jbmr.080105", "NCBI GEO series page; PubMed"],
    ["GSE13850", "smoking-related osteoporosis blood B cells", "20", "not assigned on GEO series page", "NCBI GEO Series GSE13850, Gene Expression of Circulating B Lymphocytes for Smoking-related Osteoporosis in Postmenopausal Females.", "", "NCBI GEO series page; no publication citation displayed"],
    ["GSE11917", "coronary smooth-muscle calcification", "105", "20665672", "Shalhoub V et al. Chondro/osteoblastic and cardiovascular gene modulation in human artery smooth muscle cells that calcify in the presence of phosphate and calcitriol or paricalcitol. J Cell Biochem. 2010;111:911-921.", "10.1002/jcb.22779", "NCBI GEO series page; PubMed"],
    ["GSE37558", "calcifying vascular smooth muscle cells versus osteoblasts", "32", "25380738", "Alves RD et al. Calcifying vascular smooth muscle cells and osteoblasts: independent cell types exhibiting extracellular matrix and biomineralization-related mimicries. BMC Genomics. 2014;15:965.", "10.1186/1471-2164-15-965", "NCBI GEO series page; PubMed"],
    ["GSE80614", "osteogenic/adipogenic hMSC time course", "69", "28344004", "van de Peppel J et al. Identification of Three Early Phases of Cell-Fate Determination during Osteogenic and Adipogenic Differentiation by Transcription Factor Dynamics. Stem Cell Reports. 2017;8:947-960.", "10.1016/j.stemcr.2017.02.018", "NCBI GEO series page; PubMed"],
    ["GSE159677", "calcified carotid plaque scRNA-seq", "6 samples; 3 patients", "36224302", "Alsaigh T et al. Decoding the transcriptome of calcified atherosclerotic plaque at single-cell resolution. Commun Biol. 2022;5:1084.", "10.1038/s42003-022-04056-7", "NCBI GEO series page; PubMed"],
    ["METHOD-STROBE-MR", "reporting guideline", "20-item checklist", "", "Skrivankova VW et al. Strengthening the Reporting of Observational Studies in Epidemiology Using Mendelian Randomization. JAMA. 2021.", "10.1001/jama.2021.18236", "STROBE-MR official statement"],
    ["METHOD-STROBE-MR-E&E", "reporting guideline explanation", "20-item explanation and elaboration", "", "Skrivankova VW et al. STROBE-MR: explanation and elaboration. BMJ. 2021;375:n2233.", "10.1136/bmj.n2233", "STROBE-MR official E&E"],
    ["METHOD-MR-EGGER", "MR-Egger", "method", "", "Bowden J et al. Assessing the suitability of summary data for two-sample Mendelian randomization analyses using MR-Egger regression. Int J Epidemiol. 2016;45:1961-1974.", "10.1093/ije/dyw220", "PubMed/Crossref"],
    ["METHOD-MVMR", "multivariable MR", "method", "", "Sanderson E et al. An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. Int J Epidemiol. 2019;48:713-727.", "10.1093/ije/dyy262", "OpenAlex/PubMed metadata"],
    ["METHOD-COLOC", "colocalization", "method", "", "Giambartolomei C et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet. 2014;10:e1004383.", "10.1371/journal.pgen.1004383", "Crossref/PubMed"],
    ["METHOD-PLINK", "LD clumping", "method/software", "17701901", "Purcell S et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet. 2007;81:559-575.", "10.1086/519795", "PubMed"],
    ["METHOD-1000G", "EUR reference panel", "method/resource", "26432245", "1000 Genomes Project Consortium et al. A global reference for human genetic variation. Nature. 2015;526:68-74.", "10.1038/nature15393", "PubMed"],
    ["METHOD-MSIGDB", "Hallmark gene sets", "gene-set resource", "", "Liberzon A et al. The Molecular Signatures Database Hallmark Gene Set Collection. Cell Syst. 2015;1:417-425.", "10.1016/j.cels.2015.12.004", "OpenAlex/Crossref"],
    ["METHOD-MRBASE", "summary-statistics MR platform", "method/resource", "", "Hemani G et al. The MR-Base platform supports systematic causal inference across the human phenome. eLife. 2018;7:e34408.", "10.7554/eLife.34408", "OpenAlex/Crossref"],
    ["METHOD-MEFF", "effective number of tests", "multiplicity method", "", "Li MX, Ji L. Evaluating the effective numbers of independent tests and significant p-value thresholds in commercial genotyping arrays and public imputation reference datasets. Hum Genet. 2012;131:747-756.", "10.1007/s00439-011-1118-2", "Crossref metadata"],
]


def build_reference_table():
    rows = [{"accession": r[0], "role_or_trait": r[1], "sample_size": r[2], "pmid": r[3], "citation": r[4], "doi": r[5], "verification": r[6]} for r in REFERENCE_ROWS]
    write_csv(OUT / "tables" / "Table_S1_accession_pmid_doi.csv", rows)


def build_table_manifest():
    rows = [
        {"supplement": "Table S2", "file": "tables/Table_S2_eBMD_instruments.csv", "source": "results/tables/02-mr/instruments_eBMD.csv", "description": "600 LD-independent heel-ultrasound eBMD instruments and strength metrics"},
        {"supplement": "Table S3", "file": "tables/Table_S3_finngen_endpoint_ladder.csv", "source": "results/tables/14-cardiac-ladder/ladder.csv", "description": "FinnGen cardiac ladder endpoint estimates"},
        {"supplement": "Table S4", "file": "tables/Table_S4_external_cohorts.csv", "source": "results/tables/14-cardiac-ladder/big_cohorts_all.csv", "description": "AF, HF and cardiomyopathy cohort-level estimates"},
        {"supplement": "Table S5", "file": "tables/Table_S5_circulating_trait_mvmr.csv", "source": "results/tables/06-bone-axis/axis_mvmr.csv", "description": "Ten prespecified circulating-trait MVMR models"},
        {"supplement": "Table S6", "file": "tables/Table_S6_body_composition_mvmr.csv", "source": "results/tables/11-body-composition/bodycomp_mvmr.csv", "description": "Body-composition MVMR models"},
        {"supplement": "Table S7", "file": "tables/Table_S7_colocalization.csv", "source": "results/tables/09-integration/coloc_summary.csv", "description": "Colocalization summary and ABO exclusion"},
        {"supplement": "Table S8", "file": "tables/Table_S8_pathway_geneset.csv", "source": "results/tables/12-pathway-stratified/geneset_fwer.csv; results/tables/12-pathway-stratified/pathway_contributions.csv", "description": "Pathway contribution and gene-set family-wise scan"},
        {"supplement": "Table S9", "file": "tables/Table_S9_robustness.csv", "source": "results/tables/13-robustness/negative_controls.csv; results/tables/14-cardiac-ladder/effective_tests.csv; results/tables/14-cardiac-ladder/panel_checks.csv", "description": "Negative controls, Steiger/shared-locus and effective-test records"},
        {"supplement": "Table S10", "file": "tables/Table_S10_nhanes.csv", "source": "results/tables/01-nhanes/*.csv", "description": "NHANES flow, weighted regression and mortality results"},
        {"supplement": "Table S11", "file": "tables/Table_S11_bulk_deg.csv", "source": "results/tables/03-bulk-deg/deg_summary.csv; intersection_grid.csv; pathway_convergence.csv; within_side_replication.csv", "description": "Bulk transcriptomic convergence and replication summary"},
        {"supplement": "Table S12", "file": "tables/Table_S12_scrna.csv", "source": "results/tables/07-scrna/*.csv; results/tables/09-integration/celltype_enrichment.csv", "description": "Descriptive plaque single-cell composition and programme scores"},
        {"supplement": "Table S13", "file": "tables/Table_S13_dxa_locus_matching.csv", "source": "results/tables/13-robustness/dxa_matched_loci.csv; dxa_positional_concordance.csv; dxa_concordance.csv", "description": "DXA/ultrasound locus and positional matching"},
        {"supplement": "Table S14", "file": "tables/Table_S14_pooled_af_hf.csv", "source": "results/tables/14-cardiac-ladder/big_cohorts_pooled_all_v1.csv", "description": "Durable AF/HF pooled estimates"},
    ]
    write_csv(OUT / "tables" / "Table_S0_supplement_manifest.csv", rows)


def build_csv_tables():
    copies = {
        "Table_S2_eBMD_instruments.csv": "results/tables/02-mr/instruments_eBMD.csv",
        "Table_S3_finngen_endpoint_ladder.csv": "results/tables/14-cardiac-ladder/ladder.csv",
        "Table_S4_external_cohorts.csv": "results/tables/14-cardiac-ladder/big_cohorts_all.csv",
        "Table_S5_circulating_trait_mvmr.csv": "results/tables/06-bone-axis/axis_mvmr.csv",
        "Table_S6_body_composition_mvmr.csv": "results/tables/11-body-composition/bodycomp_mvmr.csv",
        "Table_S7_colocalization.csv": "results/tables/09-integration/coloc_summary.csv",
        "Table_S8_pathway_contributions.csv": "results/tables/12-pathway-stratified/pathway_contributions.csv",
        "Table_S8_geneset_fwer.csv": "results/tables/12-pathway-stratified/geneset_fwer.csv",
        "Table_S9_negative_controls.csv": "results/tables/13-robustness/negative_controls.csv",
        "Table_S9_effective_tests.csv": "results/tables/14-cardiac-ladder/effective_tests.csv",
        "Table_S9_panel_checks.csv": "results/tables/14-cardiac-ladder/panel_checks.csv",
        "Table_S10_sample_flow.csv": "results/tables/01-nhanes/table0_sample_flow.csv",
        "Table_S10_survey_logistic.csv": "results/tables/01-nhanes/table3_survey_logistic_cvd.csv",
        "Table_S10_survey_cox.csv": "results/tables/01-nhanes/table4_survey_cox_mortality.csv",
        "Table_S10_outcome_comparison.csv": "results/tables/01-nhanes/table12_outcome_comparison.csv",
        "Table_S10_outcome_heterogeneity.csv": "results/tables/01-nhanes/table13_outcome_heterogeneity.csv",
        "Table_S11_deg_summary.csv": "results/tables/03-bulk-deg/deg_summary.csv",
        "Table_S11_intersection_grid.csv": "results/tables/03-bulk-deg/intersection_grid.csv",
        "Table_S11_pathway_convergence.csv": "results/tables/03-bulk-deg/pathway_convergence.csv",
        "Table_S11_within_side_replication.csv": "results/tables/03-bulk-deg/within_side_replication.csv",
        "Table_S12_celltype_composition.csv": "results/tables/07-scrna/celltype_composition.csv",
        "Table_S12_programme_scores.csv": "results/tables/07-scrna/programme_scores.csv",
        "Table_S12_osteo_by_tissue.csv": "results/tables/07-scrna/osteo_by_tissue.csv",
        "Table_S12_celltype_enrichment.csv": "results/tables/09-integration/celltype_enrichment.csv",
        "Table_S13_dxa_matched_loci.csv": "results/tables/13-robustness/dxa_matched_loci.csv",
        "Table_S13_dxa_positional_concordance.csv": "results/tables/13-robustness/dxa_positional_concordance.csv",
        "Table_S13_dxa_concordance.csv": "results/tables/13-robustness/dxa_concordance.csv",
        "Table_S14_pooled_af_hf.csv": "results/tables/14-cardiac-ladder/big_cohorts_pooled_all_v1.csv",
    }
    for dest, src in copies.items():
        copy_table(src, dest)


def environment_text() -> str:
    plink = ROOT / "data" / "raw" / "reference" / "plink" / "plink.exe"
    pver = subprocess.run([str(plink), "--version"], capture_output=True, text=True, check=True).stdout.strip()
    sha = hashlib.sha256(plink.read_bytes()).hexdigest()
    return f"""# Supplementary environment and reproducibility record

Generated: {datetime.now(timezone.utc).isoformat()}

## Observed software

- Python: 3.12.14 (`{sys.executable}`)
- R: 4.6.1 (recorded in `{ENV.relative_to(ROOT).as_posix()}`)
- PLINK: {pver}
- PLINK executable: `data/raw/reference/plink/plink.exe`
- PLINK SHA-256: `{sha}`
- 1000 Genomes reference: phase 3 EUR, `data/raw/reference/EUR.*`
- ReportLab: {__import__('reportlab').Version}

## Data status

- Planned downloads: 96
- Present and size-validated: 91
- Missing/unavailable in the frozen inventory: 5
- Missing/unavailable entries: 3 upstream-unavailable GSE159677 single-sample files and 2 cardiomyopathy inputs not present in the frozen inventory.
- Result-table presence is not treated as proof of complete raw-input reconstruction.
- The complete pipeline has not been independently rerun in the current project.

## Primary exposure sample-size audit

- The current harmonised eBMD file contains `n=445921` for all 12,003,351 rows.
- GWAS Catalog study page and OpenGWAS metadata report a discovery/sample_size value of `583314`.
- Difference: 137,393 people; 23.55% of the Catalog metadata value.
- Interpretation used in the manuscript: 445,921 is the analysis-file/sample-size field used for instruments and downstream calculations; 583,314 is the public study/discovery metadata value. The project does not infer that they are the same analytic cohort or a simple exclusion count.
"""


def build_supp_methods():
    text = """# Supplementary Methods

## Study design and reporting standard

This study used a layered two-sample Mendelian randomization design with an NHANES population-context analysis, a FinnGen R12 endpoint ladder, external atrial-fibrillation/heart-failure/cardiomyopathy datasets, exposure-modality triangulation, multivariable MR, colocalization, pathway/gene-set analyses, bulk transcriptomics, descriptive plaque single-cell analysis and prespecified negative controls. The reporting structure follows the 20-item STROBE-MR checklist. The checklist is supplied as a separate CSV and PDF.

## Primary exposure and the 445,921/583,314 discrepancy

The primary exposure was GCST90029004 heel-ultrasound bone density. The downloaded harmonised summary-statistics file contains a constant `n` field of 445,921 across 12,003,351 rows. The GWAS Catalog study page and OpenGWAS metadata describe the study discovery/sample size as 583,314 European-ancestry individuals. These values differ by 137,393 (23.55% of 583,314). They are retained as distinct metadata fields: 445,921 is the analysis-file/sample-size value used by the current project for instrument R虏/F-statistic calculations; 583,314 is the public study/discovery description. The records do not establish that the difference represents a simple participant exclusion, and the manuscript does not label one value as a corrected version of the other.

Six hundred genome-wide significant variants were LD-clumped using PLINK v1.90b7.2 and a 1000 Genomes phase 3 EUR reference panel (`r虏<0.001`, 10,000 kb). The retained instruments explained 危R虏=0.182, with minimum F-statistic=23.4.

## MR estimation

Outcome summary statistics were harmonised to the exposure effect allele. Palindromic variants with intermediate allele frequency were excluded. The primary estimator was random-effects IVW. MR-Egger intercept, leave-one-out, Steiger filtering, shared-locus exclusion and effective-test correction were prespecified sensitivity analyses. Multivariable models were interpreted as attenuation/residual association under the tested covariate set, not as proof that all unmeasured pathways were absent.

## Endpoint localization and pooling

The FinnGen R12 ladder used the same exposure instruments across endpoints. Eight non-anchor endpoints were evaluated with Bonferroni 伪=0.00625; an expanded 20-endpoint correlation structure yielded Meff=14.124 and 伪Meff鈮?.00354. AF and HF pooled estimates were calculated using random-effects inverse-variance pooling with DerSimonian鈥揕aird 蟿虏 and are now persisted in Table S14. The ultrasound exposure files represent multiple analyses of the UK Biobank source population and are not independent-cohort replications.

## NHANES and transcriptomic analyses

NHANES 2005鈥?018 analyses used survey weights, strata and PSUs. The final analytic sample was n=14,260 after age, femoral-neck BMD, cardiovascular-status and mortality-eligibility filters. Mortality used the NCHS linked mortality file. Bulk GEO analyses used limma with dataset-specific grouping and a paired design for GSE43292. GSE51495 was excluded after species review identified baboon samples. GSE159677 single-cell results are descriptive, based on three patients and carotid plaque rather than myocardium.

## Reproducibility boundary

The supplementary tables are derived from existing project result tables and are not a claim that all 96 raw inputs can be reconstructed in the current frozen state. Five planned inputs remain missing/unavailable, and the complete pipeline has not been independently rerun. Exact software and external-tool records are in `environment_record.txt`.
"""
    (OUT / "Supplementary_Methods.md").write_text(text, encoding="utf-8")


def build_checklist():
    rows = [
        ["1", "Name Mendelian randomization in title/abstract", "Title; Abstract", "covered", "The title and abstract name a multi-cohort Mendelian randomization study."],
        ["2", "Background and rationale", "Introduction paragraphs 1鈥?", "covered", "Observational overlap, reverse causation and shared causes are stated."],
        ["3", "Specific causal hypotheses", "Introduction paragraph 3; Supplementary Methods", "covered", "Heel-ultrasound bone density 鈫?cardiac endpoints and localization hypotheses are explicit."],
        ["4", "Study design and data sources", "Methods 2.1鈥?.3; Supplementary Table S1", "covered", "Accession, release and source mapping are supplied."],
        ["5", "Participants/samples and eligibility", "Methods 2.1鈥?.3, 2.8; Tables S1/S3/S4/S10", "covered", "GWAS sample descriptions, ancestry, endpoint counts and NHANES flow are supplied."],
        ["6", "Variables and phenotype definitions", "Methods 2.2鈥?.8; Tables S1/S3", "covered", "Exposure modality, endpoint labels, phenocodes and sample-size fields are documented."],
        ["7", "Instrument selection and strength", "Methods MR estimation; Table S2", "covered", "p threshold, LD parameters, instrument count, 危R虏 and minimum F are stated."],
        ["8", "Assumptions and pleiotropy", "Methods; Discussion; Table S9", "covered", "Relevance, independence, exclusion restriction and pleiotropy diagnostics are stated."],
        ["9", "Statistical methods", "Methods; Supplementary Methods; Tables S3鈥揝14", "covered", "IVW, MR-Egger, MVMR, colocalization, pooling and permutation parameters are recorded."],
        ["10", "Sample overlap and population structure", "Methods 2.3; Limitations; Tables S1/S4", "covered", "AF MVP independence and HF overlap boundaries are explicit."],
        ["11", "Primary estimate and uncertainty", "Results 3.1; Figure 2; Tables S4/S14", "covered", "Primary AF/HF estimates and CIs are linked to durable tables."],
        ["12", "Sensitivity analyses", "Results 3.4鈥?.7; Table S9", "covered", "Negative controls, Steiger/shared-locus and Meff outputs are supplied."],
        ["13", "Multiple testing", "Methods; Results; Tables S3/S8/S9", "covered", "Bonferroni and effective-test thresholds are stated consistently."],
        ["14", "Missing data and unavailable inputs", "Methods 2.9; Data availability; environment_record.txt", "covered with explicit boundary", "91/96 present and 5 missing/unavailable are reported; complete rerun is not claimed."],
        ["15", "Study limitations", "Discussion; Limitations", "covered", "Phenotype, overlap, ancestry, power and tissue limitations are placed with interpretation."],
        ["16", "Interpretation", "Discussion; Conclusion", "covered", "MR-consistent language is used and treatment claims are excluded."],
        ["17", "Generalisability", "Discussion; Limitations", "covered", "Ancestry, endpoint-definition, cohort-dependence and tissue boundaries are stated."],
        ["18", "Funding and conflicts", "Declarations", "author input required", "Funding and COI statements require the authors' factual entries."],
        ["19", "Data and code availability", "Declarations; Supplementary Methods", "author input required", "Public repository, commit/DOI and final access wording require the authors' entries."],
        ["20", "Reporting checklist and supplementary materials", "Supplementary files", "covered; author sign-off required", "Checklist and supplementary files are generated; authors must approve the final upload version."],
    ]
    write_csv(OUT / "STROBE-MR_checklist.csv", [{"item": r[0], "requirement": r[1], "location": r[2], "status": r[3], "evidence_or_action": r[4]} for r in rows])


def build_pdf():
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="TitleCenter", parent=styles["Title"], alignment=TA_CENTER, fontName="Helvetica-Bold", fontSize=16, leading=20, spaceAfter=10))
    styles.add(ParagraphStyle(name="H1x", parent=styles["Heading1"], fontName="Helvetica-Bold", fontSize=13, leading=16, spaceBefore=8, spaceAfter=5))
    styles.add(ParagraphStyle(name="H2x", parent=styles["Heading2"], fontName="Helvetica-Bold", fontSize=10.5, leading=13, spaceBefore=6, spaceAfter=3))
    styles.add(ParagraphStyle(name="Bodyx", parent=styles["BodyText"], fontName="Helvetica", fontSize=8.5, leading=11.5, spaceAfter=4))
    styles.add(ParagraphStyle(name="Smallx", parent=styles["BodyText"], fontName="Helvetica", fontSize=7.3, leading=9, spaceAfter=2))
    styles.add(ParagraphStyle(name="Footx", parent=styles["BodyText"], fontName="Helvetica", fontSize=7, leading=8, textColor=colors.grey))

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(colors.grey)
        canvas.drawString(18 * mm, 10 * mm, "JBMR supplementary material 鈥?evidence-backed project build")
        canvas.drawRightString(192 * mm, 10 * mm, f"Page {doc.page}")
        canvas.restoreState()

    path = OUT / "Supplementary_Methods_and_Results.pdf"
    doc = BaseDocTemplate(str(path), pagesize=A4, rightMargin=16 * mm, leftMargin=16 * mm, topMargin=15 * mm, bottomMargin=17 * mm, title="Supplementary Methods and Results")
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="normal")
    doc.addPageTemplates([PageTemplate(id="main", frames=[frame], onPage=footer)])
    story = [Paragraph("Supplementary Methods and Results", styles["TitleCenter"]), Paragraph("Heel-ultrasound bone density and non-coronary cardiac disease: a multi-cohort Mendelian randomization study", styles["Bodyx"]), Paragraph("This file is generated from the current project result tables and audit records. It is not a claim of complete raw-input reconstruction.", styles["Footx"]), Spacer(1, 5 * mm)]
    story.append(Paragraph("S1. Study design, sample-size audit and software", styles["H1x"]))
    for para in (OUT / "Supplementary_Methods.md").read_text(encoding="utf-8").split("\n\n"):
        if para.startswith("# ") or para.startswith("## "):
            continue
        if para.strip():
            story.append(Paragraph(tidy(para).replace("`", ""), styles["Bodyx"]))
    story.append(PageBreak())
    story.append(Paragraph("S2. Key pooled estimates", styles["H1x"]))
    pooled = read_csv("results/tables/14-cardiac-ladder/big_cohorts_pooled_all_v1.csv")
    data = [["Endpoint", "k", "尾", "95% CI", "P", "I虏"]]
    for r in pooled:
        data.append([r["endpoint"], r["k"], f"{float(r['beta']):.3f}", f"{float(r['lo']):.3f} to {float(r['hi']):.3f}", f"{float(r['p']):.2g}", f"{float(r['i2']):.2f}"])
    t = Table(data, colWidths=[28 * mm, 12 * mm, 20 * mm, 40 * mm, 25 * mm, 20 * mm], repeatRows=1)
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#D9EAF7")), ("GRID", (0, 0), (-1, -1), 0.25, colors.grey), ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"), ("FONTSIZE", (0, 0), (-1, -1), 8), ("ALIGN", (1, 1), (-1, -1), "RIGHT")]))
    story.append(t)
    story.append(Spacer(1, 4 * mm))
    story.append(Paragraph("The pooled AF estimate is based on five outcome datasets/definitions; the pooled HF estimate is based on four cohorts. The 'cases_total_reported' field in the source table is an audit total of reported case counts, not a deduplicated participant count.", styles["Smallx"]))
    story.append(Paragraph("S3. Evidence boundary", styles["H1x"]))
    story.append(Paragraph("The primary exposure uses the 445,921 analysis-file n field. The 583,314 value remains the public discovery/sample-size metadata value. The difference is documented, not silently harmonised. Five planned inputs remain missing or unavailable, and the complete pipeline has not been independently rerun.", styles["Bodyx"]))
    doc.build(story)
    return path


def build_strobe_pdf():
    # Use the same PDF engine for a readable upload-ready checklist.
    rows = read_csv("results/reports/JBMR-STROBE-MR浜ゅ弶鏍告煡琛╛v1.csv")
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="T", parent=styles["Title"], alignment=TA_CENTER, fontSize=15, leading=18))
    styles.add(ParagraphStyle(name="B", parent=styles["BodyText"], fontSize=7.2, leading=9))
    path = OUT / "STROBE-MR_checklist.pdf"
    doc = BaseDocTemplate(str(path), pagesize=A4, rightMargin=10 * mm, leftMargin=10 * mm, topMargin=12 * mm, bottomMargin=14 * mm, title="STROBE-MR checklist")
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="normal")
    doc.addPageTemplates([PageTemplate(id="main", frames=[frame])])
    story = [Paragraph("STROBE-MR Checklist 鈥?project crosswalk", styles["T"]), Paragraph("Official framework: 20 recommended items for transparent reporting of Mendelian randomization studies. Official site: strobe-mr.org. Current status distinguishes covered items from author-input declarations.", styles["B"]), Spacer(1, 3 * mm)]
    data = [["Item", "Requirement", "Manuscript/supplement location", "Status", "Action/evidence"]]
    for r in rows:
        data.append([r["STROBE-MR_item"], r["reporting_requirement"], r["planned_location_in_JBMR_draft"], r["current_status"], r["author_or_analysis_action"]])
    wrapped = [[Paragraph(tidy(x), styles["B"]) for x in row] for row in data]
    t = Table(wrapped, colWidths=[10 * mm, 43 * mm, 43 * mm, 25 * mm, 55 * mm], repeatRows=1)
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#D9EAF7")), ("GRID", (0, 0), (-1, -1), 0.25, colors.grey), ("VALIGN", (0, 0), (-1, -1), "TOP"), ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"), ("FONTSIZE", (0, 0), (-1, -1), 7)]))
    story.append(t)
    doc.build(story)
    return path


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    (OUT / "tables").mkdir(parents=True)
    build_reference_table()
    build_table_manifest()
    build_csv_tables()
    build_supp_methods()
    build_checklist()
    env = environment_text()
    (OUT / "environment_record.txt").write_text(env, encoding="utf-8")
    pdf = build_pdf()
    strobe_pdf = build_strobe_pdf()
    # Copy a stable audit snapshot; it is a source record, not a generated claim.
    shutil.copy2(ENV, OUT / "environment_record_project_snapshot.txt")
    manifest = []
    for f in sorted(OUT.rglob("*")):
        if f.is_file():
            manifest.append({"path": f.relative_to(OUT).as_posix(), "bytes": f.stat().st_size, "sha256": hashlib.sha256(f.read_bytes()).hexdigest()})
    (OUT / "MANIFEST_SHA256.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"created {len(manifest)} supplementary files")
    print(pdf.relative_to(ROOT).as_posix())
    print(strobe_pdf.relative_to(ROOT).as_posix())


if __name__ == "__main__":
    main()

