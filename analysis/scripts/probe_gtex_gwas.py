"""Verify GTEx v10 scope and GWAS Catalog harmonised sumstats inventory."""
import gzip
import io
import json
import urllib.error
import urllib.request

UA = {"User-Agent": "Mozilla/5.0 (research-data-probe)"}


def get(url, timeout=120, decode=True):
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = r.read()
            return r.status, (data.decode("utf-8", "replace") if decode else data)
    except urllib.error.HTTPError as e:
        return e.code, str(e)
    except Exception as e:
        return "ERR", str(e)


print("== GTEx dataset metadata ==")
st, b = get("https://gtexportal.org/api/v2/metadata/dataset")
if st == 200:
    for d in json.loads(b):
        print(f"  {d['datasetId']:10s} samples={d.get('sampleCount')} tissues={d.get('tissueCount')} "
              f"donors={d.get('subjectCount')} gencode={d.get('gencodeVersion')} dbgap={d.get('dbgapId')}")

print("\n== GTEx v10 tissue list (bulk-gex bucket XML listing) ==")
st, b = get("https://storage.googleapis.com/storage/v1/b/adult-gtex/o?prefix=bulk-gex/v10/rna-seq/&maxResults=20")
print("  ", st)
if st == 200:
    for i in json.loads(b).get("items", []):
        print("   ", i["name"], i.get("size"))
st, b = get("https://storage.googleapis.com/storage/v1/b/adult-gtex/o?delimiter=/&prefix=bulk-gex/v10/")
if st == 200:
    js = json.loads(b)
    print("  v10 prefixes:", js.get("prefixes"))

print("\n== GWAS Catalog harmonised sumstats list ==")
for u in [
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/harmonised_list.txt.gz",
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/harmonised_list.txt",
]:
    st, b = get(u)
    print("  ", st, u.split("/")[-1], (len(b) if st == 200 else str(b)[:100]))
    if st == 200:
        if u.endswith(".gz"):
            txt = gzip.decompress(b.encode("latin-1")).decode("utf-8", "replace")
        else:
            txt = b
        lines = txt.splitlines()
        print("   lines:", len(lines), "| header:", lines[0][:150])
        for kw in ("bone mineral density", "osteoporosis", "coronary artery disease", "heel"):
            hits = [l for l in lines if kw in l.lower()]
            print(f"   '{kw}': {len(hits)}")
            for h in hits[:6]:
                print("      ", h[:160])
        break

