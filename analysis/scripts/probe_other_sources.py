"""Probe availability of the non-GEO public resources needed for the
epidemiology / clinical-phenotype / shared-pathway layers.

Read-only HEAD/GET probes; prints status + a small sample of the metadata so the
report only claims what actually resolved.
"""
import json
import urllib.error
import urllib.request

UA = {"User-Agent": "Mozilla/5.0 (research-data-probe)"}


def get(url, timeout=45, binary=False):
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = r.read()
            return r.status, (len(data) if binary else data[:200000].decode("utf-8", "replace"))
    except urllib.error.HTTPError as e:
        return e.code, str(e)
    except Exception as e:
        return "ERR", str(e)


def head(url, timeout=45):
    req = urllib.request.Request(url, headers=UA, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.headers.get("Content-Length", "")
    except urllib.error.HTTPError as e:
        return e.code, str(e)
    except Exception as e:
        return "ERR", str(e)


print("=" * 70)
print("1) NHANES DXA (bone) + linked mortality (epidemiology layer)")
nhanes = {
    "2005-06 DXA femur/hip": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/DXX_D.xpt",
    "2007-08 DXA": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/DXX_E.xpt",
    "2009-10 DXA": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/DXX_F.xpt",
    "2013-14 DXA": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/DXX_H.xpt",
    "mortality 2005-06": "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHANES_2005_2006_MORT_2019_PUBLIC.dat",
    "mortality 2007-08": "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHANES_2007_2008_MORT_2019_PUBLIC.dat",
}
for k, u in nhanes.items():
    print(f"  {k:24s} {head(u)}")

print("=" * 70)
print("2) GTEx tissue catalogue (shared-pathway / normal reference layer)")
st, body = get("https://gtexportal.org/api/v2/dataset/tissueInfo")
if st == 200:
    js = json.loads(body)
    names = [d.get("tissueSiteDetailId") for d in js.get("data", [])]
    print(f"  GTEx API ok, {len(names)} tissues")
    print("  relevant:", [t for t in names if t and any(k in t for k in ("Artery", "Muscle", "Adipose", "Whole_Blood", "Cells"))])
else:
    print("  GTEx API:", st, body[:200])
st, body = head("https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz")
print("  GTEx v10 gene TPM matrix:", st, body)

print("=" * 70)
print("3) GWAS Catalog (genetic layer: BMD / osteoporosis / CAD)")
for trait in ("bone mineral density", "osteoporosis", "coronary artery disease"):
    q = urllib.parse_quote = urllib.request.quote(trait)
    st, body = get(f"https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait?diseaseTrait={q}&size=5")
    if st == 200:
        js = json.loads(body)
        print(f"  '{trait}': {js.get('page', {}).get('totalElements', '?')} studies")
        for s in js.get("_embedded", {}).get("studies", [])[:5]:
            print(f"     {s.get('accessionId')}  n={str(s.get('initialSampleSize'))[:60]}  {str(s.get('publicationInfo',{}).get('title'))[:70]}")
    else:
        print(f"  '{trait}': {st} {str(body)[:120]}")

print("=" * 70)
print("4) FinnGen R12 summary statistics (replication / MR layer)")
for prefix in ("summary_stats/M13_OSTEOPOROSIS", "summary_stats/I9_CORONARY", "summary_stats/M13_OSTEO"):
    st, body = get(f"https://storage.googleapis.com/storage/v1/b/finngen-public-data-r12/o?prefix={prefix}&maxResults=8")
    if st == 200:
        js = json.loads(body)
        items = [i["name"] for i in js.get("items", [])]
        print(f"  prefix {prefix}: {len(items)} files e.g. {items[:3]}")
    else:
        print(f"  prefix {prefix}: {st} {str(body)[:150]}")

print("=" * 70)
print("5) Open Targets Platform (shared target/evidence layer, GraphQL)")
oq = '{"query":"{ target(ensemblId: \\"ENSG00000157764\\"){ approvedSymbol } }"}'
req = urllib.request.Request("https://api.platform.opentargets.org/api/v4/graphql",
                            data=oq.encode(), headers={**UA, "Content-Type": "application/json"})
try:
    with urllib.request.urlopen(req, timeout=45) as r:
        print("  Open Targets API:", r.status, r.read()[:120].decode())
except Exception as e:
    print("  Open Targets API:", e)

print("=" * 70)
print("6) ARCHS4 / Enrichr / STRING (pathway + network layer)")
for k, u in {
    "Enrichr addList": "https://maayanlab.cloud/Enrichr/addList",
    "STRING API": "https://string-db.org/api/json/get_string_ids?identifiers=SPP1&species=9606",
    "KEGG REST": "https://rest.kegg.jp/get/hsa04380",
    "Reactome content": "https://reactome.org/ContentService/data/query/R-HSA-1474244",
}.items():
    print(f"  {k:18s} {head(u)}")

