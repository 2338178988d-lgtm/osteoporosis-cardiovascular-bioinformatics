"""GEO candidate discovery for the osteoporosis x cardiovascular comorbidity study.

Queries NCBI E-utilities (db=gds) with curated term groups, then resolves each hit
through esummary to collect accession, platform, organism, sample count and design.

Outputs: data/processed/geo_candidates/geo_hits_raw.csv
Read-only against NCBI; no bulk downloads.
"""
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path("data/processed/geo_candidates")
OUT.mkdir(parents=True, exist_ok=True)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# term groups: label -> list of alternative query strings (already GEO-syntax)
QUERIES = {
    "comorbid_direct": [
        '(osteoporosis[All Fields]) AND (cardiovascular[All Fields] OR atherosclerosis[All Fields] OR "vascular calcification"[All Fields])',
        '(osteoporosis[All Fields]) AND ("coronary"[All Fields] OR "myocardial"[All Fields] OR "heart failure"[All Fields])',
        '(bone mineral density[All Fields]) AND (atherosclerosis[All Fields] OR cardiovascular[All Fields])',
    ],
    "osteoporosis_human": [
        '"osteoporosis"[MeSH Terms] AND "Homo sapiens"[Organism] AND (expression profiling by array[DataSet Type] OR expression profiling by high throughput sequencing[DataSet Type])',
        '("bone mineral density"[All Fields] OR "osteoporotic"[All Fields] OR "low bone mass"[All Fields]) AND "Homo sapiens"[Organism] AND (expression profiling by array[DataSet Type] OR expression profiling by high throughput sequencing[DataSet Type])',
    ],
    "bone_cell_types": [
        '("osteoblast"[All Fields] OR "osteoclast"[All Fields] OR "osteocyte"[All Fields]) AND "Homo sapiens"[Organism] AND (expression profiling by array[DataSet Type] OR expression profiling by high throughput sequencing[DataSet Type])',
        '("bone marrow mesenchymal"[All Fields] OR "BMSC"[All Fields] OR "MSC osteogenic"[All Fields]) AND "Homo sapiens"[Organism]',
    ],
    "vascular_calcification": [
        '("vascular calcification"[All Fields] OR "arterial calcification"[All Fields] OR "aortic valve calcification"[All Fields] OR "coronary calcification"[All Fields]) AND (expression profiling by array[DataSet Type] OR expression profiling by high throughput sequencing[DataSet Type])',
        '("osteogenic differentiation"[All Fields] AND "vascular smooth muscle"[All Fields]) AND (expression profiling by array[DataSet Type] OR expression profiling by high throughput sequencing[DataSet Type])',
    ],
    "atherosclerosis_human": [
        '("atherosclerosis"[MeSH Terms] OR "atherosclerotic plaque"[All Fields] OR "carotid plaque"[All Fields]) AND "Homo sapiens"[Organism] AND (expression profiling by array[DataSet Type] OR expression profiling by high throughput sequencing[DataSet Type])',
        '("coronary artery disease"[All Fields] OR "coronary atherosclerosis"[All Fields]) AND "Homo sapiens"[Organism] AND expression profiling by high throughput sequencing[DataSet Type]',
    ],
    "single_cell_vascular": [
        '("single cell"[All Fields] OR "single-cell"[All Fields] OR "scRNA-seq"[All Fields]) AND (atherosclerosis[All Fields] OR "vascular"[All Fields] OR "artery"[All Fields]) AND "Homo sapiens"[Organism]',
        '("single cell"[All Fields] OR "scRNA-seq"[All Fields]) AND ("bone marrow"[All Fields] OR osteoporosis[All Fields] OR osteoblast[All Fields]) AND "Homo sapiens"[Organism]',
    ],
}


def eutil(endpoint, params):
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as exc:  # transient NCBI throttling
            if attempt == 2:
                print(f"  ! failed {endpoint}: {exc}")
                return {}
            time.sleep(2 * (attempt + 1))


def search(term, retmax=60):
    js = eutil("esearch.fcgi", {"db": "gds", "term": term, "retmax": retmax, "retmode": "json"})
    return js.get("esearchresult", {}).get("idlist", [])


def summary(uids):
    if not uids:
        return []
    js = eutil("esummary.fcgi", {"db": "gds", "id": ",".join(uids), "retmode": "json"})
    res = js.get("result", {})
    return [res[u] for u in res.get("uids", []) if u in res]


rows, seen = [], {}
for group, terms in QUERIES.items():
    for term in terms:
        uids = search(term)
        print(f"[{group}] {len(uids):3d} hits :: {term[:80]}")
        for rec in summary(uids):
            acc = rec.get("accession", "")
            if not acc.startswith("GSE"):
                continue
            if acc in seen:
                if group not in seen[acc]["groups"]:
                    seen[acc]["groups"].append(group)
                continue
            n = int(rec.get("n_samples", 0) or 0)
            row = {
                "accession": acc,
                "groups": [group],
                "title": (rec.get("title", "") or "").replace("\n", " "),
                "type": rec.get("gdstype", ""),
                "organism": rec.get("taxon", ""),
                "n_samples": n,
                "platform": ";".join(g.get("accession", "") for g in rec.get("gpl", "") or []) if isinstance(rec.get("gpl"), list) else "",
                "pubmed": rec.get("pubmedids", [""])[0] if rec.get("pubmedids") else "",
                "date": rec.get("pdat", ""),
                "summary": (rec.get("summary", "") or "").replace("\n", " ")[:1500],
            }
            if "Homo sapiens" not in row["organism"] and "single_cell" not in group:
                row["organism"] = row["organism"]  # keep non-human for cross-species validation
            seen[acc] = row
            rows.append(row)
        time.sleep(0.4)

rows.sort(key=lambda r: (-len(r["groups"]), -r["n_samples"]))
with open(OUT / "geo_hits_raw.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["accession", "groups", "title", "type", "organism", "n_samples", "platform", "pubmed", "date", "summary"])
    w.writeheader()
    for r in rows:
        r = dict(r)
        r["groups"] = "|".join(r["groups"])
        w.writerow(r)

print(f"\nunique GSE records: {len(rows)} -> {OUT/'geo_hits_raw.csv'}")
print("\nTop human GSE by group co-occurrence:")
for r in rows[:30]:
    if "Homo sapiens" in r["organism"]:
        print(f"  {r['accession']}  n={r['n_samples']:4d}  {r['type'][:28]:28s} {'|'.join(r['groups'])[:40]:40s} {r['title'][:70]}")

