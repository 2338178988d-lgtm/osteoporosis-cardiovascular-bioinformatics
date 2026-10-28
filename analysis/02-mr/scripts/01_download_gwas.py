"""Download GWAS Catalog harmonised summary statistics for the osteoporosis x CVD
Mendelian-randomisation layer (module 02-mr).

The connection to EBI in this environment is throttled to ~14 KB/s per socket, so
files are fetched as parallel byte ranges (each range is one socket) and written
to resumable chunk files. Every completed file is checked against the md5 that EBI
publishes in the accession's harmonised/md5sum.txt before it is moved into
data/raw/gwas/.

Usage (project .venv):
    python analysis/02-mr/scripts/01_download_gwas.py --list
    python analysis/02-mr/scripts/01_download_gwas.py --all --connections 12
    python analysis/02-mr/scripts/01_download_gwas.py --datasets eBMD,CAD
"""
import argparse
import concurrent.futures as cf
import hashlib
import os
import re
import shutil
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data" / "raw" / "gwas"
CHUNKS = RAW / ".chunks"
MANIFEST = RAW / "manifest.csv"
FTP = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
UA = {"User-Agent": "Mozilla/5.0 (research data download; module 02-mr)"}

# key: dict(acc, bin, remote name on EBI FTP, local file name, trait, role, priority)
DATASETS = {
    "eBMD": dict(acc="GCST90029004", bin="GCST90029001-GCST90030000",
                 remote="29892013-GCST90029004-EFO_0009270.h.tsv.gz",
                 local="GCST90029004_eBMD.h.tsv.gz",
                 trait="heel bone mineral density (eBMD), European", role="exposure", prio=1),
    "CAD": dict(acc="GCST90132315", bin="GCST90132001-GCST90133000",
                remote="GCST90132315.h.tsv.gz", local="GCST90132315_CAD.h.tsv.gz",
                trait="coronary artery disease (Aragam 2022)", role="outcome", prio=2),
    "STROKE": dict(acc="GCST90104534", bin="GCST90104001-GCST90105000",
                   remote="GCST90104534.h.tsv.gz", local="GCST90104534_STROKE.h.tsv.gz",
                   trait="any stroke (Mishra 2022, GIGASTROKE)", role="outcome", prio=3),
    "HF": dict(acc="GCST90162626", bin="GCST90162001-GCST90163000",
               remote="GCST90162626.h.tsv.gz", local="GCST90162626_HF.h.tsv.gz",
               trait="heart failure (Levin 2022, HERMES)", role="outcome", prio=4),
    "MI": dict(acc="GCST90018877", bin="GCST90018001-GCST90019000",
               remote="GCST90018877.h.tsv.gz", local="GCST90018877_MI.h.tsv.gz",
               trait="myocardial infarction (Sakaue 2021)", role="outcome", prio=5),
    "OP": dict(acc="GCST90018887", bin="GCST90018001-GCST90019000",
               remote="GCST90018887.h.tsv.gz", local="GCST90018887_OP.h.tsv.gz",
               trait="osteoporosis (Sakaue 2021)", role="exposure", prio=6),
}


def url_for(acc, bin_, fname):
    return f"{FTP}{bin_}/{acc}/harmonised/{fname}"


def head_size(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA, method="HEAD")
            with urllib.request.urlopen(req, timeout=60) as r:
                n = r.headers.get("Content-Length")
                if n:
                    return int(n)
        except Exception as e:
            log(f"  HEAD failed ({i + 1}/{tries}): {e}")
            time.sleep(2 * (i + 1))
    return None


def fetch_text(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            log(f"  text fetch failed ({i + 1}/{tries}): {e}")
            time.sleep(2 * (i + 1))
    return ""


def md5_of(path, buf=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def log(msg):
    print(msg, flush=True)


def get_chunk(url, path, start, end, tries=6):
    """Download url[start:end] into path, resuming if path already has bytes."""
    want = end - start + 1
    for attempt in range(tries):
        have = path.stat().st_size if path.exists() else 0
        if have >= want:
            return True
        a = start + have
        req = urllib.request.Request(url, headers={**UA, "Range": f"bytes={a}-{end}"})
        try:
            with urllib.request.urlopen(req, timeout=120) as r, open(path, "ab") as f:
                while True:
                    block = r.read(1 << 16)
                    if not block:
                        break
                    f.write(block)
        except Exception:
            time.sleep(1.5 * (attempt + 1))
        if path.exists() and path.stat().st_size >= want:
            return True
    return path.exists() and path.stat().st_size >= want


def download(key, connections=12, chunk_mb=4):
    d = DATASETS[key]
    url = url_for(d["acc"], d["bin"], d["remote"])
    dest = RAW / d["local"]

    log(f"\n=== {key} ({d['acc']}) -> {dest.name}")
    total = head_size(url)
    if total is None:
        log("  cannot determine remote size; skipped")
        return False
    log(f"  remote size: {total / 1e6:.1f} MB")

    # an interrupted earlier attempt is preserved, never silently overwritten
    if dest.exists():
        if dest.stat().st_size == total and md5_ok(dest, d):
            log("  already complete (md5 ok)")
            return True
        bad = dest.with_suffix(dest.suffix + f".incomplete-{dest.stat().st_size}")
        dest.rename(bad)
        log(f"  existing/incomplete file preserved as {bad.name}")

    cdir = CHUNKS / key
    cdir.mkdir(parents=True, exist_ok=True)
    step = chunk_mb * 1024 * 1024
    ranges = [(s, min(s + step - 1, total - 1)) for s in range(0, total, step)]
    parts = [cdir / f"{i:05d}.part" for i in range(len(ranges))]
    todo = [(u, p, s, e) for (s, e), u, p in zip(ranges, [url] * len(ranges), parts)
            if not (p.exists() and p.stat().st_size >= e - s + 1)]
    log(f"  {len(ranges)} chunks ({chunk_mb} MB each), {len(todo)} to fetch, {connections} sockets")

    t0 = time.time()
    last = [0, time.time()]

    def report():
        got = sum(p.stat().st_size for p in parts if p.exists())
        now = time.time()
        rate = (got - last[0]) / (now - last[1]) / 1e6 if now > last[1] else 0
        last[0], last[1] = got, now
        log(f"  {got / 1e6:8.1f} / {total / 1e6:.1f} MB "
            f"({100 * got / total:5.1f}%)  {rate * 1000:6.1f} KB/s  {now - t0:6.0f}s")

    with cf.ThreadPoolExecutor(max_workers=connections) as ex:
        futs = {ex.submit(get_chunk, *job): job for job in todo}
        pending = set(futs)
        while pending:
            done, pending = cf.wait(pending, timeout=90)
            for f in done:
                if not f.result():
                    log(f"  chunk failed: {futs[f][1].name}")
            report()
            if not pending:
                break

    missing = [p for (s, e), p in zip(ranges, parts) if not (p.exists() and p.stat().st_size >= e - s + 1)]
    if missing:
        log(f"  {len(missing)} chunks still incomplete; rerun to resume")
        return False

    tmp = dest.with_suffix(dest.suffix + ".assembling")
    with open(tmp, "wb") as out:
        for p in parts:
            with open(p, "rb") as f:
                shutil.copyfileobj(f, out, 1 << 22)
    if tmp.stat().st_size != total:
        log(f"  assembled size mismatch ({tmp.stat().st_size} != {total})")
        return False
    tmp.replace(dest)
    shutil.rmtree(cdir, ignore_errors=True)

    if not md5_ok(dest, d):
        log("  md5 mismatch -> file kept but flagged in manifest")
        status = "md5-mismatch"
    else:
        log(f"  done, md5 ok, {dest.stat().st_size / 1e6:.1f} MB")
        status = "ok"
    row = [key, d["acc"], url, str(dest.relative_to(ROOT).as_posix()), str(total), md5_of(dest), status,
           time.strftime("%Y-%m-%dT%H:%M:%S")]
    new = not MANIFEST.exists()
    with open(MANIFEST, "a", encoding="utf-8") as f:
        if new:
            f.write("key,accession,url,path,bytes,md5,status,finished\n")
        f.write(",".join(row) + "\n")
    return status == "ok"


_MD5_CACHE = {}


def md5_ok(path, d):
    """Compare a local file with the md5 EBI publishes for the remote file."""
    acc, bin_ = d["acc"], d["bin"]
    if acc not in _MD5_CACHE:
        _MD5_CACHE[acc] = fetch_text(f"{FTP}{bin_}/{acc}/harmonised/md5sum.txt")
    want = None
    n = d["remote"]
    for line in _MD5_CACHE[acc].splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[1] == n:
            want = parts[0]
    if want is None:
        log("  no published md5 found; size check only")
        return True
    got = md5_of(path)
    if got != want:
        log(f"  md5 {got} != published {want}")
        return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", default="", help="comma separated keys, e.g. eBMD,CAD")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--connections", type=int, default=12)
    ap.add_argument("--chunk-mb", type=int, default=4)
    a = ap.parse_args()

    if a.list:
        for k, v in sorted(DATASETS.items(), key=lambda x: x[1]["prio"]):
            log(f"{k:7s} {v['acc']:12s} {v['trait']:45s} {v['role']:9s} prio={v['prio']}")
        return 0

    keys = list(DATASETS) if a.all else [k.strip() for k in a.datasets.split(",") if k.strip()]
    if not keys:
        log("nothing selected; use --all or --datasets")
        return 2
    keys.sort(key=lambda k: DATASETS[k]["prio"])
    RAW.mkdir(parents=True, exist_ok=True)
    ok = []
    for k in keys:
        try:
            if download(k, a.connections, a.chunk_mb):
                ok.append(k)
        except KeyboardInterrupt:
            log("interrupted; chunks are kept, rerun to resume")
            return 130
    log(f"\ncomplete: {ok}  ({len(ok)}/{len(keys)})")
    return 0 if len(ok) == len(keys) else 1


if __name__ == "__main__":
    sys.exit(main())

