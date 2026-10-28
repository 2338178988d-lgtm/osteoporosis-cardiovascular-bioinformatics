"""Supplementary GEO search to fill remaining gaps:
human calcified artery tissue (ex vivo), human bone tissue transcriptome,
human bone-marrow single-cell, and CAVD valve tissue cohorts.
"""
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path("data/processed/geo_candidates")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

QUERIES = {
    "artery_calcified_tissue": '("calcification"[All Fields] OR "calcified"[All Fields]) AND ("artery"[All Fields] OR "aorta"[All Fields] OR "aortic"[All Fields]) AND "Homo sapiens"[Organism] AND (expression profiling by array[DataSet Type] OR expression profiling by high throughput sequencing[DataSet Type])',
    "bone_tissue_human": '("trabecular bone"[All Fields] OR "cortical bone"[All Fields] OR "bone biopsy"[All Fields] OR "femoral head"[All Fields]) AND "Homo sapiens"[Organism] AND (expression profiling by array[DataSet Type] OR expression profiling by high throughput sequencing[DataSet Type])',
    "bone_marrow_sc": '("bone marrow"[All Fields]) AND ("single cell"[All Fields] OR "scRNA"[All Fields] OR "single-cell"[All Fields]) AND "Homo sapiens"[Organism]',
    "cavd_valve": '("aortic valve"[All Fields] AND ("stenosis"[All Fields] OR "calcific"[All Fields])) AND "Homo sapiens"[Organism] AND (expression profiling by array[DataSet Type] OR expression profiling by high throughput sequencing[DataSet Type])',
    "osteogenic_vsmc": '("vascular smooth muscle"[All Fields] AND ("osteogenic"[All Fields] OR "osteoblast"[All Fields] OR "mineralization"[All Fields] OR "calcification"[All Fields])) AND "Homo sapiens"[Organism]',
}


def eutil(endpoint, params):
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    for i in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception:
            time.sleep(2 * (i + 1))
    return {}


rows = []
for label, term in QUERIES.items():
    js = eutil("esearch.fcgi", {"db": "gds", "term": term, "retmax": 80, "retmode": "json"})
    uids = js.get("esearchresult", {}).get("idlist", [])
    print(f"[{label}] {len(uids)} hits")
    if not uids:
        continue
    sj = eutil("esummary.fcgi", {"db": "gds", "id": ",".join(uids), "retmode": "json"})
    res = sj.get("result", {})
    for u in res.get("uids", []):
        rec = res[u]
        acc = rec.get("accession", "")
        if not acc.startswith("GSE") or "Homo sapiens" not in (rec.get("taxon") or ""):
            continue
        if int(rec.get("n_samples") or 0) < 12:
            continue
        rows.append({
            "gap": label,
            "accession": acc,
            "n_samples": rec.get("n_samples"),
            "type": rec.get("gdstype", ""),
            "date": rec.get("pdat", ""),
            "pubmed": (rec.get("pubmedids") or [""])[0],
            "title": (rec.get("title") or "").replace("\n", " ")[:160],
            "summary": (rec.get("summary") or "").replace("\n", " ")[:500],
        })
    time.sleep(0.5)

seen, uniq = set(), []
for r in sorted(rows, key=lambda r: -int(r["n_samples"])):
    if r["accession"] in seen:
        continue
    seen.add(r["accession"])
    uniq.append(r)

with open(OUT / "geo_gap_candidates.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["gap", "accession", "n_samples", "type", "date", "pubmed", "title", "summary"])
    w.writeheader()
    w.writerows(uniq)
print(f"\n{len(uniq)} unique GSE (n>=12, human) -> {OUT/'geo_gap_candidates.csv'}")
for r in uniq[:40]:
    print(f"  {r['accession']:11s} n={int(r['n_samples']):>4} [{r['gap'][:22]:22s}] {r['title'][:85]}")

