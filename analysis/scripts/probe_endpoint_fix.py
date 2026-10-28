"""Resolve the correct endpoints/paths for sources that failed the first probe."""
import json
import urllib.error
import urllib.request

UA = {"User-Agent": "Mozilla/5.0 (research-data-probe)"}


def get(url, timeout=45):
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.read()[:120000].decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        return e.code, str(e)
    except Exception as e:
        return "ERR", str(e)


print("== NHANES DXA filename variants ==")
for u in [
    "https://wwwn.cdc.gov/Nchs/Nhanes/2007-2008/DXX_E.XPT",
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/DXX_E.XPT",
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/DXX_F.XPT",
    "https://wwwn.cdc.gov/Nchs/Nhanes/2009-2010/DXX_F.XPT",
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/DXX_D.XPT",
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DXX_J.XPT",
]:
    print(" ", u.split("Nhanes/")[-1], get(u)[0])

print("\n== GTEx API tissue endpoints ==")
for u in [
    "https://gtexportal.org/api/v2/dataset/tissueInfo",
    "https://gtexportal.org/api/v2/reference/tissueInfo",
    "https://gtexportal.org/api/v2/dataset/datasetInfo",
    "https://gtexportal.org/api/v2/expression/medianGeneExpression?gencodeId=ENSG00000164692.18",
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/",
]:
    st, b = get(u)
    print(" ", st, u)
    if st == 200 and len(b) < 5000:
        print("    ", b[:400].replace("\n", " "))

print("\n== GWAS Catalog API variants ==")
for u in [
    "https://www.ebi.ac.uk/gwas/rest/api/studies?size=2",
    "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/search/findByTrait?trait=osteoporosis&size=2",
    "https://www.ebi.ac.uk/gwas/api/v2/studies?size=2",
    "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/search/findByGene?geneName=SPP1&size=2",
]:
    st, b = get(u)
    print(" ", st, u)
    if st == 200:
        print("    ", b[:300].replace("\n", " "))

print("\n== FinnGen R12 bucket layout ==")
for prefix in ["summary_stats/", "summary_stats/r12/", "", "summary_stats/finngen_R12_M13"]:
    st, b = get(f"https://storage.googleapis.com/storage/v1/b/finngen-public-data-r12/o?prefix={prefix}&maxResults=6")
    print("  prefix=", repr(prefix), st)
    if st == 200:
        js = json.loads(b)
        print("    ", [i["name"] for i in js.get("items", [])])

print("\n== FinnGen endpoint metadata (R12) ==")
st, b = get("https://r12.finngen.fi/api/endpoints")
print("  r12 endpoints:", st, b[:200].replace("\n", " "))
st, b = get("https://r12.risteys.finngen.fi/api/endpoints")
print("  risteys endpoints:", st, b[:300].replace("\n", " "))

