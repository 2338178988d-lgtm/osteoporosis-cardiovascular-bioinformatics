"""Build the JBMR reference-closure table and extend Supplementary Table S1.

All entries are based on already verified project records or live official
metadata gathered for this project. No citation is inferred from an accession
without a source record.
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
S1 = ROOT / "results/supplementary/tables/Table_S1_accession_pmid_doi.csv"
AUDIT = ROOT / "results/reports/JBMR_accession_reference_audit_v2.csv"
OUT = ROOT / "results/reports/JBMR-reference-closure_v1.csv"

ADDITIONAL_S1 = [
    {
        "accession": "GCST90000066",
        "role_or_trait": "systolic blood pressure; MVMR/diagnostic context",
        "sample_size": "810,865 European; additional multi-ancestry participants",
        "pmid": "33230300",
        "citation": "Surendran P et al. Discovery of rare variants associated with blood pressure regulation through meta-analysis of 1.3 million individuals. Nat Genet. 2020;52:1314-1332.",
        "doi": "10.1038/s41588-020-00713-x",
        "verification": "GWAS Catalog study page; PubMed metadata",
    },
    {
        "accession": "GCST005413",
        "role_or_trait": "type 2 diabetes; MVMR",
        "sample_size": "up to 12,931 cases and 57,196 controls",
        "pmid": "29358691",
        "citation": "Bon脿s-Guarch S et al. Re-analysis of public genetic data reveals a rare X-chromosomal variant associated with type 2 diabetes. Nat Commun. 2018;9:321.",
        "doi": "10.1038/s41467-017-02380-9",
        "verification": "GWAS Catalog study page; PubMed metadata",
    },
    {
        "accession": "DATA-NHANES-2005-2018",
        "role_or_trait": "NHANES public-use survey files",
        "sample_size": "2005-2018; seven two-year cycles",
        "pmid": "",
        "citation": "National Center for Health Statistics. National Health and Nutrition Examination Survey public-use data files.",
        "doi": "",
        "verification": "Official CDC/NCHS data portal: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/",
    },
    {
        "accession": "DATA-NCHS-LMF-2019",
        "role_or_trait": "NCHS 2019 public-use linked mortality file",
        "sample_size": "Linked to NHANES 1999-2018 public-use records",
        "pmid": "",
        "citation": "National Center for Health Statistics. Public-use linked mortality files, 2019; data dictionary and release documentation.",
        "doi": "",
        "verification": "Project data dictionary plus official NCHS linked-mortality directory: https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/",
    },
    {
        "accession": "DATA-FINNGEN-R12",
        "role_or_trait": "FinnGen R12 summary statistics and endpoint manifest",
        "sample_size": "Release R12; endpoint-specific",
        "pmid": "",
        "citation": "FinnGen. FinnGen R12 public summary statistics and endpoint manifest.",
        "doi": "",
        "verification": "Official manifest: https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv; accessed 2026-10-23/24",
    },
]

CLOSURE = [
    ("REF-CONDURACHE-2024", "Condurache DG et al. Bone health, cardiovascular disease, and imaging outcomes in UK Biobank: a causal analysis. JBMR Plus. 2024;8(6):ziae058.", "38784722", "10.1093/jbmrpl/ziae058", "Discussion comparator", "verified_full_text"),
    ("REF-SALIH-2025", "Salih AM et al. Bone mineral density and cardiovascular diseases: a two-sample Mendelian randomization study. JBMR Plus. 2025;9(5):ziaf037.", "40191156", "10.1093/jbmrpl/ziaf037", "Discussion comparator and power accounting", "verified_full_text"),
    ("REF-QIAN-2026", "Qian Y et al. Proxied therapeutic inhibition on sclerostin and atrial fibrillation risk. Communications Medicine. 2026;6(1).", "42218280", "10.1038/s43856-026-01694-4", "Discussion treatment-proxy boundary", "verified_full_text"),
    ("GCST90029004", "Loh PR et al. Mixed-model association for biobank-scale datasets. Nat Genet. 2018;50:906-908.", "29892013", "10.1038/s41588-018-0144-6", "Primary heel-ultrasound eBMD source", "verified_metadata"),
    ("GCST90025982", "Barton AR et al. Whole-exome imputation within UK Biobank powers rare coding variant association and fine-mapping analyses. Nat Genet. 2021;53:1260-1269.", "34226706", "10.1038/s41588-021-00892-1", "Heel-ultrasound BMD T-score source", "verified_metadata"),
    ("GCST006979", "Morris JA et al. An atlas of genetic influences on osteoporosis in humans and mice. Nat Genet. 2019;51:258-266.", "30598549", "10.1038/s41588-018-0302-x", "Heel-ultrasound BMD source", "verified_metadata"),
    ("GCST90446641/GCST90446628", "Miao J et al. Valid inference for machine-learning-assisted genome-wide association studies. Nat Genet. 2024;56:2361-2369.", "39349818", "10.1038/s41588-024-01934-0", "DXA femoral-neck BMD source", "verified_metadata"),
    ("GCST90568458", "Qian Y et al. Genome-wide association studies of over 30,000 samples with bone mineral density at multiple skeletal sites and its clinical relevance. Genomics Proteomics Bioinformatics. 2025;23(5).", "41206123", "10.1093/gpbjnl/qzaf097", "DXA femoral-neck BMD source", "verified_metadata"),
    ("GCST90162626", "Levin MG et al. Genome-wide association and multi-trait analyses characterize the common genetic architecture of heart failure. Nat Commun. 2022;13:6914.", "36376295", "10.1038/s41467-022-34216-6", "HERMES heart-failure source", "verified_metadata"),
    ("GCST009541", "Shah S et al. Genome-wide association and Mendelian randomisation analysis provide insights into the pathogenesis of heart failure. Nat Commun. 2020;11:163.", "31919418", "10.1038/s41467-019-13690-5", "Heart-failure replication source", "verified_metadata"),
    ("GCST90018806/GCST90018877/GCST90018887/GCST90018586/GCST90018667", "Sakaue S et al. A cross-population atlas of genetic associations for 220 human phenotypes. Nat Genet. 2021;53:1415-1424.", "34594039", "10.1038/s41588-021-00931-x", "Cross-population HF/MI/osteoporosis/BBJ sources", "verified_metadata"),
    ("GCST90132315", "Aragam KG et al. Discovery and systematic characterization of risk variants and genes for coronary artery disease in over a million participants. Nat Genet. 2022;54:1803-1815.", "36474045", "10.1038/s41588-022-01233-6", "CAD source", "verified_metadata"),
    ("GCST90104534", "Mishra A et al. Stroke genetics informs drug discovery and risk prediction across ancestries. Nature. 2022;611:115-123.", "36180795", "10.1038/s41586-022-05165-3", "Stroke source", "verified_metadata"),
    ("GCST90244711", "Khurshid S et al. Clinical and genetic associations of deep learning-derived cardiac magnetic resonance-based left ventricular mass. Nat Commun. 2023;14:1558.", "36944631", "10.1038/s41467-023-37173-w", "Left-ventricular-mass source", "verified_metadata"),
    ("GCST90435260", "Tadros R et al. Large-scale genome-wide association analyses identify novel genetic loci and mechanisms in hypertrophic cardiomyopathy. Nat Genet. 2025;57:530-538.", "39966646", "10.1038/s41588-025-02087-4", "Indexed LV-mass source", "verified_metadata"),
    ("GCST90000025", "Pei YF et al. The genetic architecture of appendicular lean mass characterized by association analysis in the UK Biobank study. Commun Biol. 2020;3:608.", "33097823", "10.1038/s42003-020-01334-0", "Appendicular lean mass source", "verified_metadata"),
    ("GCST90000066", "Surendran P et al. Discovery of rare variants associated with blood pressure regulation through meta-analysis of 1.3 million individuals. Nat Genet. 2020;52:1314-1332.", "33230300", "10.1038/s41588-020-00713-x", "Systolic blood pressure MVMR source", "verified_metadata"),
    ("GCST005413", "Bon脿s-Guarch S et al. Re-analysis of public genetic data reveals a rare X-chromosomal variant associated with type 2 diabetes. Nat Commun. 2018;9:321.", "29358691", "10.1038/s41467-017-02380-9", "Type 2 diabetes MVMR source", "verified_metadata"),
    ("GCST90469223/GCST90470711/GCST90468440/GCST90468997/GCST90470372", "Loya H et al. A scalable variational inference approach for increased mixed-model association power. Nat Genet. 2025;57:461-468.", "39789286", "10.1038/s41588-024-02044-7", "FGF23/SOST/BGLAP/DKK1/PTH mediator sources", "verified_metadata"),
    ("DATA-NHANES-2005-2018", "National Center for Health Statistics. National Health and Nutrition Examination Survey public-use data files.", "", "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", "NHANES population layer", "official_data_source"),
    ("DATA-NCHS-LMF-2019", "National Center for Health Statistics. Public-use linked mortality files, 2019.", "", "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/", "NHANES mortality linkage", "official_data_source"),
    ("DATA-FINNGEN-R12", "FinnGen. FinnGen R12 public summary statistics and endpoint manifest.", "", "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv", "FinnGen R12 outcome layer", "official_data_source"),
    ("METHOD-STROBE-MR", "Skrivankova VW et al. Strengthening the Reporting of Observational Studies in Epidemiology Using Mendelian Randomization. JAMA. 2021.", "", "10.1001/jama.2021.18236", "Reporting guideline", "verified_metadata"),
    ("METHOD-STROBE-MR-E&E", "Skrivankova VW et al. STROBE-MR: explanation and elaboration. BMJ. 2021;375:n2233.", "", "10.1136/bmj.n2233", "Reporting guideline explanation", "verified_metadata"),
    ("METHOD-MR-EGGER", "Bowden J et al. Assessing the suitability of summary data for two-sample Mendelian randomization analyses using MR-Egger regression. Int J Epidemiol. 2016;45:1961-1974.", "", "10.1093/ije/dyw220", "MR-Egger method", "verified_metadata"),
    ("METHOD-MVMR", "Sanderson E et al. An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. Int J Epidemiol. 2019;48:713-727.", "", "10.1093/ije/dyy262", "MVMR method", "verified_metadata"),
    ("METHOD-COLOC", "Giambartolomei C et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet. 2014;10:e1004383.", "", "10.1371/journal.pgen.1004383", "Colocalization method", "verified_metadata"),
    ("METHOD-PLINK", "Purcell S et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet. 2007;81:559-575.", "17701901", "10.1086/519795", "LD clumping software/method", "verified_metadata"),
    ("METHOD-1000G", "1000 Genomes Project Consortium et al. A global reference for human genetic variation. Nature. 2015;526:68-74.", "26432245", "10.1038/nature15393", "EUR reference panel", "verified_metadata"),
    ("METHOD-MSIGDB", "Liberzon A et al. The Molecular Signatures Database Hallmark Gene Set Collection. Cell Syst. 2015;1:417-425.", "", "10.1016/j.cels.2015.12.004", "Gene-set resource", "verified_metadata"),
    ("METHOD-MRBASE", "Hemani G et al. The MR-Base platform supports systematic causal inference across the human phenome. eLife. 2018;7:e34408.", "", "10.7554/eLife.34408", "Summary-statistics MR resource", "verified_metadata"),
    ("METHOD-MEFF", "Li MX, Ji L. Evaluating the effective numbers of independent tests and significant p-value thresholds in commercial genotyping arrays and public imputation reference datasets. Hum Genet. 2012;131:747-756.", "", "10.1007/s00439-011-1118-2", "Multiplicity method", "verified_metadata"),
]


def read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        return list(reader), list(reader.fieldnames or [])


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    rows, fields = read_csv(S1)
    existing = {row["accession"] for row in rows}
    for row in ADDITIONAL_S1:
        if row["accession"] not in existing:
            rows.append(row)
    write_csv(S1, rows, fields)

    closure_fields = ["ref_key", "citation", "pmid", "doi_or_official_url", "role", "status"]
    closure_rows = [dict(zip(closure_fields, item)) for item in CLOSURE]
    write_csv(OUT, closure_rows, closure_fields)

    print(f"updated {S1.relative_to(ROOT)} rows={len(rows)}")
    print(f"wrote {OUT.relative_to(ROOT)} rows={len(closure_rows)}")


if __name__ == "__main__":
    main()

