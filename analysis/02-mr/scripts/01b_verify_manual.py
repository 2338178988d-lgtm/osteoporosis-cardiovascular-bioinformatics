"""Check GWAS summary-statistic files that were downloaded outside this project
(browser, download manager, wget on another machine) and register them.

The files may be dropped into data/raw/gwas/ under either the original EBI file name
or the short local name used by the analysis scripts. Each file is checked for exact
byte size, the md5 that EBI publishes (where published), and a full gzip read to the
end of stream (CRC check). Only files that pass are renamed to the local name and
recorded in data/raw/gwas/manifest.csv.

Nothing is ever deleted: a file that fails stays where it is (use --replace to move a
bad copy aside into data/raw/gwas/_replaced/ when the good copy arrives).

Usage:
    python analysis/02-mr/scripts/01b_verify_manual.py --list
    python analysis/02-mr/scripts/01b_verify_manual.py
    python analysis/02-mr/scripts/01b_verify_manual.py --replace
"""
import argparse
import gzip
import hashlib
import importlib.util
import json
import shutil
import sys
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("dl", HERE / "01_download_gwas.py")
dl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dl)

ROOT = dl.ROOT
RAW = dl.RAW
MANIFEST = dl.MANIFEST
INDEX = ROOT / "data" / "processed" / "mr" / "gwas_file_index.json"
REPLACED = RAW / "_replaced"


def build_index(refresh=False):
    if INDEX.exists() and not refresh:
        return json.loads(INDEX.read_text(encoding="utf-8"))
    rows = []
    for key, d in sorted(dl.DATASETS.items(), key=lambda x: x[1]["prio"]):
        md5 = None
        txt = dl.fetch_text(f"{dl.FTP}{d['bin']}/{d['acc']}/harmonised/md5sum.txt")
        for line in txt.splitlines():
            p = line.split()
            if len(p) == 2 and p[1] == d["remote"]:
                md5 = p[0]
        rows.append(dict(key=key, accession=d["acc"], local=d["local"], remote=d["remote"],
                         bytes=dl.head_size(dl.url_for(d["acc"], d["bin"], d["remote"])),
                         md5=md5, url=dl.url_for(d["acc"], d["bin"], d["remote"])))
    INDEX.parent.mkdir(parents=True, exist_ok=True)
    INDEX.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    return rows


def md5_of(path, buf=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def gzip_ok(path):
    """Read the whole stream: catches truncation (EOF before end of stream) and CRC errors."""
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 22):
                pass
        return True, ""
    except Exception as e:
        return False, f"{type(e).__name__}: {e}"


def find_candidates(rec):
    remote = rec.get("remote") or dl.DATASETS[rec["key"]]["remote"]
    names = {rec["local"], remote}
    return [p for n in names for p in [RAW / n] if p.exists()]


def check(rec, replace=False):
    key, want = rec["key"], rec["bytes"]
    cands = find_candidates(rec)
    if not cands:
        print(f"{key:7s} MISSING   expected {rec['local']} ({want / 1e6:.1f} MB)")
        return False
    for path in cands:
        size = path.stat().st_size
        if size != want:
            print(f"{key:7s} BAD SIZE  {path.name}: {size / 1e6:.1f} MB, expected {want / 1e6:.1f} MB "
                  f"({want - size:+d} bytes) -> incomplete")
            continue
        got = md5_of(path)
        if rec["md5"] and got != rec["md5"]:
            print(f"{key:7s} BAD MD5   {path.name}: {got}, published {rec['md5']} -> corrupt")
            continue
        ok, msg = gzip_ok(path)
        if not ok:
            print(f"{key:7s} BAD GZIP  {path.name}: {msg}")
            continue
        if path.name != rec["local"]:
            target = RAW / rec["local"]
            if target.exists():
                if not replace:
                    print(f"{key:7s} CONFLICT  {target.name} already exists; rerun with --replace")
                    continue
                REPLACED.mkdir(parents=True, exist_ok=True)
                stale = REPLACED / f"{target.name}.{int(time.time())}"
                target.replace(stale)
                print(f"           stale {target.name} moved to {stale.relative_to(ROOT).as_posix()}")
            path.replace(target)
            print(f"{key:7s} renamed   {path.name} -> {target.name}")
            path = target
        note = "md5+size+gzip ok" if rec["md5"] else "size+gzip ok (EBI publishes no md5)"
        print(f"{key:7s} OK        {path.name}  {size / 1e6:.1f} MB  ({note})")
        write_manifest(rec, path)
        return True
    return False


def write_manifest(rec, path):
    line_key = rec["accession"]
    rows = []
    header = "key,accession,url,path,bytes,md5,status,finished\n"
    if MANIFEST.exists():
        rows = [r for r in MANIFEST.read_text(encoding="utf-8").splitlines()[1:] if r.strip()]
    rows = [r for r in rows if not r.split(",")[1:2] == [line_key]]
    row = ",".join([rec["key"], rec["accession"], rec["url"],
                    path.relative_to(ROOT).as_posix(), str(path.stat().st_size),
                    md5_of(path), "ok-manual" if rec["md5"] else "ok-manual-nomd5",
                    time.strftime("%Y-%m-%dT%H:%M:%S")])
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(header + "\n".join(rows + [row]) + "\n", encoding="utf-8")


def show_list(recs):
    print("Put each file in  data/raw/gwas/  under either name below, then run this script\n")
    for r in recs:
        print(f"[{r['key']}] {r['local']}   (or {r['remote']})")
        print(f"    size {r['bytes'] / 1e6:.1f} MB   md5 {r['md5'] or 'not published by EBI'}")
        print(f"    {r['url']}")
    print(f"\ntotal {sum(r['bytes'] for r in recs) / 1e9:.2f} GB")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true", help="print URLs, names, sizes and md5")
    ap.add_argument("--refresh", action="store_true", help="re-query sizes and md5 from EBI")
    ap.add_argument("--replace", action="store_true",
                    help="move a stale same-name file aside before placing the verified copy")
    a = ap.parse_args()

    recs = build_index(a.refresh)
    if a.list:
        show_list(recs)
        return 0
    if not RAW.exists():
        print(f"{RAW} does not exist")
        return 2
    ok = [r["key"] for r in recs if check(r, a.replace)]
    missing = [r["key"] for r in recs if r["key"] not in ok]
    print(f"\nverified: {len(ok)}/{len(recs)}")
    if missing:
        print(f"still needed: {', '.join(missing)}")
    return 0 if not missing else 1


if __name__ == "__main__":
    sys.exit(main())

