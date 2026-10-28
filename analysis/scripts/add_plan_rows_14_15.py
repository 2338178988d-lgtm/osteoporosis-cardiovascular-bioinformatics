"""Append the module 14 / 15 inputs to results/tables/download_plan.csv (idempotent).

module 14 (cardiac ladder) fetches eight FinnGen R12 endpoints, but never materialises the
~815 MB summary-statistics file for each: FinnGen publishes a tabix index (`.tbi`, ~1.7 MB)
beside every file, so the ~600 instrument positions are read over HTTP range requests
(`analysis/scripts/_tabix_http.py`) and cached locally. The plan therefore carries the index,
not the full file.

module 15 (modality triangulation) needs full files, because the instrument sets for the new
exposures have to be built by scanning the whole file for p < 5e-8: two further DXA femoral
neck files and two heel-ultrasound files from publications other than the one behind the
current exposure, so that "modality" and "source" can be told apart.

Byte sizes are read from the server at write time so the downloader's size check is exact.
"""
import csv
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PLAN = ROOT / "results" / "tables" / "download_plan.csv"
FG = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/"
GC = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
UA = {"User-Agent": "Mozilla/5.0 (research data download)"}


def gwas_url(g, tail):
    n = int(g[4:])
    lo = (n // 1000) * 1000
    return f"{GC}GCST{lo + 1:06d}-GCST{lo + 1000:06d}/{g}/{tail}"


def size(url):
    """Size via a ranged GET: the GWAS FTP answers 404 to HEAD but serves the file normally."""
    req = urllib.request.Request(url, headers={**UA, "Range": "bytes=0-0"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return int(r.headers["Content-Range"].rsplit("/", 1)[1])


# (key, phenocode, cases, note)
FG_ADD = [
    ("FG_AF", "I9_AF", 63_532, "心房颤动（FinnGen R12）：63,532 例 — 非冠脉、非心室的心律通路"),
    ("FG_HYPERTENS", "I9_HYPTENS", 154_630, "高血压（FinnGen R12）：154,630 例 — 负荷本身，阶梯底端"),
    ("FG_HYPERTENS_HD", "I9_HYPTENSHD", 11_728, "高血压性心脏病（FinnGen R12）：11,728 例 — 负荷与心肌之间"),
    ("FG_NONISCH_CM", "I9_NONISCHCARDMYOP", 13_917, "非缺血性心肌病（FinnGen R12）：13,917 例 — 心衰无冠心病的可用替身"),
    ("FG_CARDMYO", "I9_CARDMYO", 7_649, "心肌病（FinnGen R12）：7,649 例 — 敏感性"),
    ("FG_CAVS", "I9_CAVS_OPERATED", 12_418, "钙化性主动脉瓣狭窄手术（FinnGen R12）：12,418 例 — 钙化悖论"),
    ("FG_VHD", "I9_VHD", 28_602, "瓣膜性心脏病（FinnGen R12）：28,602 例 — 结构性心脏病总入口"),
    ("FG_MI_STRICT", "I9_MI_STRICT", 31_666, "心肌梗死（严格，FinnGen R12）：31,666 例 — 冠脉对照"),
]

GWAS_ADD = [
    ("FNBMD_2", "GCST90446628", "harmonised/GCST90446628.h.tsv.gz",
     "DXA 股骨颈 BMD 第二份（GCST90446628, 47,205；与 GCST90446641 同一篇 2024 方法学论文，技术重复）"),
    ("FNBMD_3", "GCST90568458", "harmonised/GCST90568458.h.tsv.gz",
     "DXA 股骨颈 BMD 第三份（GCST90568458, 32,017；PMID 41206123，独立于 2024 方法学论文）"),
    ("QUS_2", "GCST90025982", "harmonised/34226706-GCST90025982-EFO_0009270.h.tsv.gz",
     "足跟超声 BMD T 值（GCST90025982, 445,855；PMID 34226706，独立于当前暴露的 Morris 文件）"),
    ("QUS_3", "GCST006979", "harmonised/30598549-GCST006979-EFO_0009270.h.tsv.gz",
     "足跟超声 BMD（GCST006979, 426,824 英国人；PMID 30598549，第三个独立来源）"),
]


def main():
    plan = pd.read_csv(PLAN)
    have = set(plan["dataset_id"])
    rows = []
    for key, pheno, cases, note in FG_ADD:
        if key in have:
            continue
        url = f"{FG}finngen_R12_{pheno}.gz.tbi"
        rows.append(dict(dataset_id=key, category="ladder", url=url,
                         local_path=f"data/raw/finngen/finngen_R12_{pheno}.gz.tbi",
                         bytes=size(url), md5="-",
                         note=f"{note}；tabix 索引，按位点远程取数"))
    for key, acc, tail, note in GWAS_ADD:
        if key in have:
            continue
        url = gwas_url(acc, tail)
        rows.append(dict(dataset_id=key, category="robustness", url=url,
                         local_path=f"data/raw/gwas/{key}_{acc}.h.tsv.gz",
                         bytes=size(url), md5="-", note=note))
    if not rows:
        print("plan already carries every row; nothing to add")
        return 0
    out = pd.concat([plan, pd.DataFrame(rows)], ignore_index=True)
    out.to_csv(PLAN, index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"added {len(rows)} rows; plan now has {len(out)}")
    for r in rows:
        print(f"  {r['dataset_id']:16s} {r['category']:10s} {r['bytes'] / 1e6:8.1f} MB  {r['local_path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
