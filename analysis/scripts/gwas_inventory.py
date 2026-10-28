"""Inventory GWAS Catalog studies for the traits this comorbidity study needs,
keeping genome-wide-scale studies and flagging summary-statistics availability.

Outputs: data/processed/gwas_candidates/gwas_studies.csv
"""
import csv
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path("data/processed/gwas_candidates")
OUT.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0 (research-data-probe)"}
BASE = "https://www.ebi.ac.uk/gwas/api/v2/studies"

TRAITS = [
    "bone mineral density", "heel bone mineral density", "osteoporosis",
    "coronary artery disease", "coronary heart disease", "myocardial infarction",
    "vascular calcification", "aortic valve stenosis", "carotid intima media thickness",
    "arterial stiffness", "pulse wave velocity", "osteoporotic fracture",
]


def get_json(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
                return json.loads(r.read().decode("utf-8", "replace"))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return {}
            time.sleep(2)
        except Exception:
            time.sleep(2)
    return {}


def parse_n(desc_list):
    """Largest N parsed out of GWAS Catalog free-text sample descriptions."""
    best = 0
    for d in desc_list or []:
        for m in re.findall(r"([\d,]{4,})\s*(?:individuals|subjects|participants|cases|samples|people|women|men|adults)", d):
            try:
                best = max(best, int(m.replace(",", "")))
            except ValueError:
                pass
        for m in re.findall(r"([\d,]{4,})\s*(?:European|East Asian|African|Hispanic|South Asian)", d):
            try:
                best = max(best, int(m.replace(",", "")))
            except ValueError:
                pass
    return best


rows = []
for trait in TRAITS:
    page, total, collected = 0, None, 0
    while True:
        url = f"{BASE}?efoTrait={urllib.parse.quote(trait)}&size=100&page={page}"
        js = get_json(url)
        studies = js.get("_embedded", {}).get("studies", [])
        if total is None:
            total = js.get("page", {}).get("totalElements", 0)
        if not studies:
            break
        for s in studies:
            n = parse_n(s.get("initialSampleDescription"))
            rows.append({
                "queried_trait": trait,
                "accession": s.get("accessionId"),
                "reported_trait": s.get("reportedTrait"),
                "efo_terms": "|".join(f"{t['key']}:{t['label']}" for t in s.get("efoTraits", [])),
                "n_discovery": n,
                "ancestry": "|".join(s.get("discoverySampleAncestry") or []),
                "summary_stats": s.get("summaryStatistics"),
                "assoc_count": s.get("associationCount"),
                "pubmed": s.get("pubmedId"),
                "first_author": s.get("firstAuthor"),
                "year": (s.get("publicationDate") or "")[:4],
                "journal": s.get("journal"),
                "title": (s.get("title") or "")[:220],
            })
        collected += len(studies)
        page += 1
        time.sleep(0.3)
        if collected >= min(total, 400):
            break
    print(f"{trait:34s} total={total:5d} collected={collected:4d}")

# dedupe by accession (a study can appear under several queried traits)
seen, uniq = set(), []
for r in rows:
    if r["accession"] in seen:
        continue
    seen.add(r["accession"])
    uniq.append(r)

keys = ["queried_trait", "accession", "reported_trait", "efo_terms", "n_discovery", "ancestry",
        "summary_stats", "assoc_count", "pubmed", "first_author", "year", "journal", "title"]
with open(OUT / "gwas_studies.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    w.writerows(uniq)

print(f"\nunique studies: {len(uniq)} -> {OUT/'gwas_studies.csv'}")
big = [r for r in uniq if r["n_discovery"] >= 20000 and r["summary_stats"] not in (None, "", "NA")]
big.sort(key=lambda r: -r["n_discovery"])
print(f"\n>=20k samples AND summary stats available: {len(big)}")
for r in big[:25]:
    print(f"  {r['accession']:12s} n={r['n_discovery']:>7} {r['reported_trait'][:52]:52s} {r['first_author'][:16]:16s} {r['year']}")

