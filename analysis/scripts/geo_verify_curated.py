"""Verify curated GEO candidates: pull real metadata (title, n_samples, platform,
organism, PMID, design summary) for a hand-picked accession list.

Outputs: data/processed/geo_candidates/geo_curated_verified.csv
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

CURATED = {
    "bone_bulk": [
        "GSE56815", "GSE56814", "GSE56816", "GSE7158", "GSE13850", "GSE7429",
        "GSE35956", "GSE35958", "GSE51495", "GSE80614", "GSE157322",
    ],
    "plaque_bulk": [
        "GSE100927", "GSE43292", "GSE28829", "GSE57691", "GSE23561",
        "GSE221911", "GSE310095", "GSE218474",
    ],
    "calcification": [
        "GSE37558", "GSE84012", "GSE148219", "GSE285688", "GSE288136",
        "GSE284683", "GSE309462",
    ],
    "single_cell": [
        "GSE159677", "GSE131780", "GSE155512", "GSE224273", "GSE131778",
    ],
    "extra_candidates": [
        "GSE317068", "GSE159121", "GSE271848", "GSE199718", "GSE247957",
        "GSE120521", "GSE10394", "GSE20129",
    ],
}


def eutil(endpoint, params):
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as exc:
            if attempt == 2:
                print(f"  ! {endpoint} failed: {exc}")
                return {}
            time.sleep(2 * (attempt + 1))


def fetch(acc):
    js = eutil("esearch.fcgi", {"db": "gds", "term": f"{acc}[Accession]", "retmode": "json"})
    ids = js.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None
    js2 = eutil("esummary.fcgi", {"db": "gds", "id": ids[0], "retmode": "json"})
    res = js2.get("result", {})
    for u in res.get("uids", []):
        rec = res[u]
        if rec.get("accession") == acc:
            return rec
    return None


rows = []
for cat, accs in CURATED.items():
    for acc in accs:
        rec = fetch(acc)
        time.sleep(0.4)
        if not rec:
            print(f"{acc:12s} [{cat}] NOT FOUND")
            rows.append({"category": cat, "accession": acc, "status": "not_found"})
            continue
        gpls = rec.get("gpl", "")
        if isinstance(gpls, list):
            gpls = ";".join(g.get("accession", "") for g in gpls)
        row = {
            "category": cat,
            "accession": acc,
            "status": "ok",
            "title": (rec.get("title") or "").replace("\n", " "),
            "n_samples": rec.get("n_samples", ""),
            "type": rec.get("gdstype", ""),
            "organism": rec.get("taxon", ""),
            "platform": gpls,
            "pubmed": (rec.get("pubmedids") or [""])[0],
            "date": rec.get("pdat", ""),
            "summary": (rec.get("summary") or "").replace("\n", " ")[:900],
        }
        rows.append(row)
        print(f"{acc:12s} [{cat}] n={row['n_samples']:>4} {row['organism'][:14]:14s} {row['type'][:26]:26s} {row['title'][:78]}")

keys = ["category", "accession", "status", "title", "n_samples", "type", "organism", "platform", "pubmed", "date", "summary"]
with open(OUT / "geo_curated_verified.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)
print(f"\nwrote {len(rows)} rows -> {OUT/'geo_curated_verified.csv'}")

