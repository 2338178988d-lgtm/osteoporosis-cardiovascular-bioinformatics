"""Prepare the large-cohort replication inputs (plan v2, section 4).

Adds the AF_US file to the download plan as a full download (579 MB, no tabix index
published) and pulls just the tabix indexes for the three indexed files.

Why not download the others in full: GCST90624411 is 1.7 GB, GCST90473571 and GCST90473575 are
4.0 GB each. Their harmonised releases carry .tbi indexes, so reading 600 positions over HTTP
range requests costs a few hundred KB instead of 13 GB.
"""
import csv
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PLAN = ROOT / "results" / "tables" / "download_plan.csv"
BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
UA = {"User-Agent": "Mozilla/5.0 (research data download)"}
IDX = ROOT / "data" / "raw" / "gwas" / "_tbi"

AF_US = ("AF_US", "GCST90475970", "GCST90475001-GCST90476000", "GCST90475970.tsv.gz",
         "房颤，美国欧洲裔（MVP PheCode 427.21，71,577 例 / 366,195 对照，PMID 39024449）"
         "— 非芬兰非 UKB，独立性最好")

INDEXED = [
    ("AF_BIG", "GCST90624411", "GCST90624001-GCST90625000",
     "harmonised/GCST90624411.h.tsv.gz.tbi",
     "房颤，欧洲 252,438 例（PMID 40645996）— 最大样本，但含芬兰与英国"),
    ("AF_NFE", "GCST90473571", "GCST90473001-GCST90474000",
     "harmonised/GCST90473571.h.tsv.gz.tbi",
     "房颤，非芬兰欧洲 40,273 例（UKB ICD10 I48，PMID 40770095）— 与芬盖零重叠"),
    ("HF_BIG", "GCST90455657", "GCST90455001-GCST90456000",
     "harmonised/GCST90455657.h.tsv.gz.tbi",
     "心衰，欧洲 168,363 例（PMID 40195560）— 比 HERMES 大 3.5 倍"),
    ("HF_NFE", "GCST90473575", "GCST90473001-GCST90474000",
     "harmonised/GCST90473575.h.tsv.gz.tbi",
     "心衰，非芬兰欧洲 19,629 例（UKB ICD10 I50，PMID 40770095）— 与芬盖零重叠"),
]


def ranged_size(url):
    req = urllib.request.Request(url, headers={**UA, "Range": "bytes=0-0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        cr = r.headers.get("Content-Range")
        return int(cr.rsplit("/", 1)[1]) if cr and "/" in cr else int(r.headers["Content-Length"])


def main():
    plan = pd.read_csv(PLAN)
    have = set(plan["dataset_id"])
    rows = []
    key, acc, bucket, main_name, note = AF_US
    if key not in have:
        url = f"{BASE}{bucket}/{acc}/{main_name}"
        rows.append(dict(dataset_id=key, category="ladder", url=url,
                         local_path=f"data/raw/gwas/{key}_{acc}.tsv.gz",
                         bytes=ranged_size(url), md5="-", note=note))
    for key, acc, bucket, tbi, note in INDEXED:
        if key in have:
            continue
        url = f"{BASE}{bucket}/{acc}/{tbi}"
        rows.append(dict(dataset_id=key, category="ladder", url=url,
                         local_path=f"data/raw/gwas/_tbi/{acc}.tbi",
                         bytes=ranged_size(url), md5="-", note=note))
    if rows:
        out = pd.concat([plan, pd.DataFrame(rows)], ignore_index=True)
        out.to_csv(PLAN, index=False, quoting=csv.QUOTE_MINIMAL)
        print(f"added {len(rows)} rows; plan now has {len(out)}")
        for r in rows:
            print(f"  {r['dataset_id']:8s} {r['bytes'] / 1e6:9.1f} MB  {r['local_path']}")
    else:
        print("plan already carries every row")

    # the tabix indexes are small enough to fetch inline
    IDX.mkdir(parents=True, exist_ok=True)
    for _, acc, bucket, tbi, _ in INDEXED:
        dest = IDX / f"{acc}.tbi"
        if dest.exists():
            print(f"  index {acc}.tbi present")
            continue
        data = urllib.request.urlopen(
            urllib.request.Request(f"{BASE}{bucket}/{acc}/{tbi}", headers=UA),
            timeout=180).read()
        dest.write_bytes(data)
        print(f"  index {acc}.tbi {len(data) / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
