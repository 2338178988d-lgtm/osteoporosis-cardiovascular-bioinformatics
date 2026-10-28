"""Build probe -> gene symbol maps for the GEO platforms used by 03-bulk-deg.

Two sources are used, in this order:
  1. the platform annotation file `annot/GPLxxx.annot.gz` on the GEO FTP (small, already
     has a "Gene symbol" column);
  2. the platform table from GEO's text interface (`acc.cgi?...&view=data`) for platforms
     that publish no annotation file; there the symbol comes from `GENE_SYMBOL` or, for
     transcript-array designs, from the second field of `gene_assignment`.

Output: data/processed/geo/_platform/<GSE>_symbols.tsv with columns probe, symbol, one
file per series (several series share a platform, so the platform table is cached once
in data/raw/geo/_platform/).

Usage:
    python analysis/03-bulk-deg/scripts/00_build_symbols.py
    python analysis/03-bulk-deg/scripts/00_build_symbols.py --gse GSE56815,GSE100927
"""
import argparse
import gzip
import io
import re
import sys
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
GEO = ROOT / "data" / "processed" / "geo"
RAW_PLAT = ROOT / "data" / "raw" / "geo" / "_platform"
OUT_PLAT = GEO / "_platform"
UA = {"User-Agent": "Mozilla/5.0 (research; module 03-bulk-deg)"}


def ftp_stub(acc):
    n = int(acc[3:])
    return f"GPL{n // 1000}nnn" if n >= 1000 else "GPLnnn"


def fetch(url, dest, timeout=300):
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
        data = r.read()
    dest.write_bytes(data)
    return dest


def from_annot(acc):
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{ftp_stub(acc)}/{acc}/annot/{acc}.annot.gz"
    dest = RAW_PLAT / f"{acc}.annot.gz"
    try:
        fetch(url, dest)
    except Exception:
        return None
    with gzip.open(dest, "rt", encoding="utf-8", errors="replace") as f:
        # the file opens with ^ / ! / # comment lines; the table header follows
        head = None
        for line in f:
            if line[:1] in ("^", "!", "#"):
                continue
            head = line.rstrip("\n").split("\t")
            break
        if not head:
            return None
        idx = {c: i for i, c in enumerate(head)}
        pcol = next((c for c in ("ID", "Probe Set ID", "SPOT_ID") if c in idx), head[0])
        scol = next((c for c in ("Gene symbol", "GENE_SYMBOL", "Gene Symbol") if c in idx), None)
        if scol is None:
            return None
        rows = []
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(idx[pcol], idx[scol]):
                continue
            rows.append((p[idx[pcol]], p[idx[scol]]))
    d = pd.DataFrame(rows, columns=["probe", "symbol"])
    return d[d["symbol"].str.strip() != ""]


def from_platform_table(acc):
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=data"
    dest = RAW_PLAT / f"{acc}_platform_table.txt"
    try:
        fetch(url, dest, timeout=600)
    except Exception:
        return None
    text = dest.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"!platform_table_begin\n(.*?)\n!platform_table_end", text, re.S)
    if not m:
        return None
    tbl = pd.read_csv(io.StringIO(m.group(1)), sep="\t", dtype=str, low_memory=False)
    pcol = "ID"
    if "GENE_SYMBOL" in tbl:
        sym = tbl["GENE_SYMBOL"].fillna("")
    elif "gene_assignment" in tbl:
        # "NM_001 // SYMBOL // description // cytoband // entrez /// ..."
        sym = tbl["gene_assignment"].fillna("").map(
            lambda s: s.split(" // ")[1] if " // " in s else "")
    else:
        return None
    d = pd.DataFrame({"probe": tbl[pcol], "symbol": sym.str.strip()})
    d = d[(d["symbol"] != "") & (d["symbol"] != "---")]
    return d.drop_duplicates("probe")


def build(acc):
    for fn in (from_annot, from_platform_table):
        d = fn(acc)
        if d is not None and len(d):
            return d.drop_duplicates("probe")
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gse", default="")
    a = ap.parse_args()
    OUT_PLAT.mkdir(parents=True, exist_ok=True)
    series = ([g.strip() for g in a.gse.split(",") if g.strip()] if a.gse
              else sorted(p.name for p in GEO.iterdir() if p.is_dir() and p.name.startswith("GSE")))
    cache, summary = {}, []
    for g in series:
        pf = GEO / g / "phenotype.csv"
        if not pf.exists():
            print(f"{g}: not processed yet, skipped")
            continue
        plats = pd.read_csv(pf)["platform"].dropna().unique()
        if not len(plats):
            print(f"{g}: no platform recorded, skipped")
            continue
        acc = plats[0]
        if acc not in cache:
            print(f"{acc}: fetching platform annotation ...", flush=True)
            cache[acc] = build(acc)
        d = cache[acc]
        if d is None:
            print(f"{g} ({acc}): no symbol source available, probes stay as they are")
            continue
        out = OUT_PLAT / f"{g}_symbols.tsv"
        d.to_csv(out, sep="\t", index=False)
        summary.append(dict(gse=g, platform=acc, probes_mapped=len(d),
                            unique_symbols=d["symbol"].nunique()))
        print(f"{g} ({acc}): {len(d)} probes -> {d['symbol'].nunique()} symbols")
    if summary:
        pd.DataFrame(summary).to_csv(OUT_PLAT / "symbols_summary.csv", index=False)


if __name__ == "__main__":
    sys.exit(main())

