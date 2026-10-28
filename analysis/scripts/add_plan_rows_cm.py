"""Add the cardiomyopathy replication inputs to the download plan.

Three files, chosen after checking what actually exists (see plan v3 section 1):

  CM_MVP2  GCST90475949  20,680 European cardiomyopathy cases (MVP PheCode), tabix index
  CM_MVP1  GCST90475951  19,708 European cardiomyopathy cases (MVP PheCode), no index
  CM_NFE   GCST90473562   3,578 non-Finnish European cardiomyopathy cases (UKB ICD),
                          tabix index - the only cardiomyopathy file with zero FinnGen overlap

GCST007714 (Gonzalez 2019, 1,816 non-ischaemic cardiomyopathy cases) is not included: at 1,816
cases it has about 15 % power for the effect of interest, which cannot decide anything.
"""
import csv
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PLAN = ROOT / "results" / "tables" / "download_plan.csv"
BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
UA = {"User-Agent": "Mozilla/5.0 (research data download)"}

ADD = [
    ("CM_MVP2", "GCST90475949", "GCST90475001-GCST90476000",
     "harmonised/GCST90475949.h.tsv.gz.tbi",
     "蹇冭倢鐥咃紝娆ф床 20,680 渚嬶紙MVP PheCode锛孭MID 39024449锛夆€?tabix 绱㈠紩锛屾寜浣嶇偣鍙栨暟"),
    ("CM_MVP1", "GCST90475951", "GCST90475001-GCST90476000",
     "GCST90475951.tsv.gz",
     "蹇冭倢鐥咃紝娆ф床 19,708 渚嬶紙MVP PheCode锛孭MID 39024449锛夆€?鏃犵储寮曪紝鏁翠唤涓嬭浇"),
    ("CM_NFE", "GCST90473562", "GCST90473001-GCST90474000",
     "harmonised/GCST90473562.h.tsv.gz.tbi",
     "蹇冭倢鐥咃紝闈炶姮鍏版娲?3,578 渚嬶紙UKB ICD锛孭MID 40770095锛夆€?涓庤姮鐩栭浂閲嶅彔锛屽姛鏁堢害 0.25"),
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
    for key, acc, bucket, tail, note in ADD:
        if key in have:
            continue
        url = f"{BASE}{bucket}/{acc}/{tail}"
        suffix = ".tbi" if tail.endswith(".tbi") else ".h.tsv.gz" if "harmonised" in tail else ".tsv.gz"
        folder = "_tbi" if suffix == ".tbi" else ""
        local = f"data/raw/gwas/{folder + '/' if folder else ''}{key}_{acc}{suffix}"
        rows.append(dict(dataset_id=key, category="ladder", url=url, local_path=local,
                         bytes=ranged_size(url), md5="-", note=note))
    if not rows:
        print("plan already carries every row")
        return 0
    out = pd.concat([plan, pd.DataFrame(rows)], ignore_index=True)
    out.to_csv(PLAN, index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"added {len(rows)} rows; plan now has {len(out)}")
    for r in rows:
        print(f"  {r['dataset_id']:8s} {r['bytes'] / 1e6:9.1f} MB  {r['local_path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

