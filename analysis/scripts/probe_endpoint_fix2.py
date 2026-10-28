"""Finish resolving: FinnGen endpoint ids, GWAS Catalog v2 trait search + sumstats,
GTEx tissue catalogue endpoint, and cohort/EHR resource availability."""
import csv
import gzip
import io
import json
import urllib.error
import urllib.request

UA = {"User-Agent": "Mozilla/5.0 (research-data-probe)"}


def get(url, timeout=90, decode=True):
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = r.read()
            return r.status, (data.decode("utf-8", "replace") if decode else data)
    except urllib.error.HTTPError as e:
        return e.code, str(e)
    except Exception as e:
        return "ERR", str(e)


print("== FinnGen R12 endpoint list: osteoporosis / CVD related ==")
st, body = get("https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv")
print("  manifest status:", st)
if st == 200:
    lines = body.splitlines()
    print("  header:", lines[0][:200])
    hits = [l for l in lines if any(k in l.upper() for k in ("OSTEOPOR", "OSTEOpen".upper(), "FRACTURE", "BONE_MIN", "M13_"))]
    print(f"  osteoporosis/bone endpoints: {len(hits)}")
    for l in hits[:15]:
        print("   ", l[:180])
    cv = [l for l in lines if any(k in l.upper() for k in ("I9_CORONARY", "I9_MI", "I9_AF", "I9_HEARTFAIL", "CORONARY_ATHEROSCLER"))]
    print(f"  CVD endpoints: {len(cv)}")
    for l in cv[:10]:
        print("   ", l[:180])

print("\n== GWAS Catalog v2: trait search params ==")
for u in [
    "https://www.ebi.ac.uk/gwas/api/v2/studies?efoTraitId=EFO_0003882&size=3",
    "https://www.ebi.ac.uk/gwas/api/v2/studies?trait=osteoporosis&size=3",
    "https://www.ebi.ac.uk/gwas/api/v2/efoTraits/EFO_0003882",
    "https://www.ebi.ac.uk/gwas/api/v2/studies/GCST006979",
]:
    st, b = get(u)
    print(" ", st, u)
    if st == 200:
        print("    ", b[:220].replace("\n", " "))

print("\n== GWAS Catalog harmonised sumstats FTP layout ==")
for u in [
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006979/",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006979/harmonised/",
]:
    st, b = get(u)
    print(" ", st, u.split("summary_statistics/")[-1])
    if st == 200:
        print("    ", b[:300].replace("\n", " "))

print("\n== GTEx tissue catalogue alternatives ==")
for u in [
    "https://gtexportal.org/api/v2/dataset/tissueInfo?datasetId=gtex_v10",
    "https://gtexportal.org/api/v2/metadata/dataset",
    "https://gtexportal.org/api/v2/expression/medianTranscriptExpression?gencodeId=ENSG00000164692.18",
    "https://gtexportal.org/api/v2/association/singleTissueEqtl?gencodeId=ENSG00000164692.18&tissueSiteDetailId=Artery_Tibial",
]:
    st, b = get(u)
    print(" ", st, u[:95])
    if st == 200:
        print("    ", b[:260].replace("\n", " "))

print("\n== Cohort / EHR resources (epidemiology + clinical phenotype) ==")
for k, u in {
    "UK Biobank showcase": "https://biobank.ndph.ox.ac.uk/showcase/",
    "dbGaP": "https://www.ncbi.nlm.nih.gov/gap/",
    "CHARLS": "http://charls.pku.edu.cn/",
    "MIMIC-IV physionet": "https://physionet.org/content/mimiciv/3.1/",
    "NHANES DXA 2013-14 codebook": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/DXX_H.htm",
    "UKB DXA field 78 (heel BMD)": "https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=78",
    "OpenGWAS": "https://api.opengwas.io/api/",
}.items():
    st, _ = get(u, timeout=45)
    print(f"  {k:30s} {st}")

