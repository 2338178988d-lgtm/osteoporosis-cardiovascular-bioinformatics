"""Download the files listed in results/tables/download_plan.csv that are not present yet.

Parallel byte ranges over several sockets (each range is its own request, which is what
gets throughput on throttled hosts), resumable: finished ranges are skipped, so rerunning
continues instead of restarting. Every file is checked for the planned byte size and, for
gzip members, read to the end of the stream before it is moved into place.

Usage:
    python analysis/scripts/download_plan.py --list
    python analysis/scripts/download_plan.py --category geo
    python analysis/scripts/download_plan.py --only GSE100927,GSE56815 --connections 8
"""
import argparse
import concurrent.futures as cf
import csv
import gzip
import shutil
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PLAN = ROOT / "results" / "tables" / "download_plan.csv"
CHUNKS = ROOT / "data" / "raw" / ".downloads"
MANIFEST = ROOT / "data" / "raw" / "download_manifest.csv"
UA = {"User-Agent": "Mozilla/5.0 (research data download)"}


def log(msg):
    print(msg, flush=True)


def head_size(url, tries=3):
    """Remote size, via HEAD and, where the host refuses HEAD, a one-byte ranged GET.

    ftp.ebi.ac.uk answers 404 to HEAD on some files while serving the same URL normally, so
    a HEAD-only probe reports "cannot determine remote size" for files that are present.
    """
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA, method="HEAD")
            with urllib.request.urlopen(req, timeout=60) as r:
                n = r.headers.get("Content-Length")
                if n:
                    return int(n)
        except Exception as e:
            log(f"  HEAD failed ({i + 1}/{tries}): {e}")
        try:
            req = urllib.request.Request(url, headers={**UA, "Range": "bytes=0-0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                cr = r.headers.get("Content-Range")          # "bytes 0-0/12345"
                if cr and "/" in cr:
                    return int(cr.rsplit("/", 1)[1])
        except Exception as e:
            log(f"  ranged GET size probe failed ({i + 1}/{tries}): {e}")
        time.sleep(2 * (i + 1))
    return None


def gzip_ok(path):
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 22):
                pass
        return True
    except Exception:
        return False


def get_range(url, path, start, end, tries=5):
    """Fetch one byte range, appending only what is still missing.

    A retry after a dropped connection appends from start+have. If the server ignores the
    Range header it answers 200 with the whole file, and appending that would silently
    corrupt the part, so a non-206 answer is discarded instead of written.
    """
    want = end - start + 1
    for attempt in range(tries):
        have = path.stat().st_size if path.exists() else 0
        if have >= want:
            return True
        req = urllib.request.Request(url, headers={**UA, "Range": f"bytes={start + have}-{end}"})
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                if r.status != 206:
                    time.sleep(1.5 * (attempt + 1))
                    continue
                with open(path, "ab") as f:
                    left = want - have
                    while left > 0:
                        block = r.read(min(1 << 16, left))
                        if not block:
                            break
                        f.write(block)
                        left -= len(block)
        except Exception:
            time.sleep(1.5 * (attempt + 1))
    return path.exists() and path.stat().st_size >= want


def download(row, connections=8, chunk_mb=2):
    url, dest = row["url"], ROOT / row["local_path"]
    planned = int(row["bytes"]) if row["bytes"] else None
    log(f"\n=== {row['dataset_id']} -> {row['local_path']}")
    if dest.exists() and planned and dest.stat().st_size == planned:
        log("  already present with the planned size; skipped")
        return True
    total = head_size(url)
    if total is None:
        log("  cannot determine remote size; skipped")
        return False
    if planned and total != planned:
        log(f"  note: remote size {total} differs from plan {planned}; using remote size")
    dest.parent.mkdir(parents=True, exist_ok=True)

    # chunks are keyed by chunk size as well as dataset: mixing runs with different
    # --chunk-mb would otherwise reuse part 00000 of one size for a range of another and
    # assemble a file that is silently the wrong length
    cdir = CHUNKS / f"{row['dataset_id']}_{chunk_mb}mb"
    cdir.mkdir(parents=True, exist_ok=True)
    step = chunk_mb * 1024 * 1024
    ranges = [(s, min(s + step - 1, total - 1)) for s in range(0, total, step)]
    parts = [cdir / f"{i:05d}.part" for i in range(len(ranges))]
    todo = [(url, p, s, e) for (s, e), p in zip(ranges, parts)
            if not (p.exists() and p.stat().st_size >= e - s + 1)]
    log(f"  {total / 1e6:.1f} MB, {len(ranges)} chunks, {len(todo)} to fetch on {connections} sockets")
    t0, last = time.time(), [0, time.time()]
    if todo:
        with cf.ThreadPoolExecutor(max_workers=connections) as ex:
            futs = [ex.submit(get_range, *job) for job in todo]
            done = 0
            for f in cf.as_completed(futs):
                done += 1
                if done % 5 == 0 or done == len(futs):
                    got = sum(p.stat().st_size for p in parts if p.exists())
                    now = time.time()
                    rate = (got - last[0]) / (now - last[1]) / 1e6 if now > last[1] else 0
                    last[0], last[1] = got, now
                    log(f"  {got / 1e6:7.1f} / {total / 1e6:.1f} MB ({100 * got / total:5.1f}%) "
                        f"{rate * 1000:7.1f} KB/s {now - t0:5.0f}s")

    missing = [p for (s, e), p in zip(ranges, parts) if not (p.exists() and p.stat().st_size >= e - s + 1)]
    if missing:
        log(f"  {len(missing)} chunks incomplete; rerun to resume")
        return False
    tmp = dest.with_suffix(dest.suffix + ".assembling")
    with open(tmp, "wb") as out:
        for p in parts:
            with open(p, "rb") as f:
                shutil.copyfileobj(f, out, 1 << 22)
    if tmp.stat().st_size != total:
        log(f"  assembled size {tmp.stat().st_size} != {total}")
        return False
    if dest.suffix in (".gz", ".tgz") and not gzip_ok(tmp):
        tmp.replace(dest.with_suffix(dest.suffix + ".corrupt"))
        log("  gzip stream is not readable end-to-end; kept as .corrupt")
        return False
    tmp.replace(dest)
    shutil.rmtree(cdir, ignore_errors=True)
    log(f"  done, {dest.stat().st_size / 1e6:.1f} MB, {time.time() - t0:.0f}s")
    write_manifest(row)
    return True


def write_manifest(row):
    dest = ROOT / row["local_path"]
    new = not MANIFEST.exists()
    rows = []
    if not new:
        with open(MANIFEST, encoding="utf-8") as f:
            rows = [r for r in list(csv.reader(f))[1:] if r and r[0] != row["dataset_id"]]
    with open(MANIFEST, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["dataset_id", "category", "local_path", "bytes", "md5", "verified_at"])
        for r in rows:
            w.writerow(r)
        w.writerow([row["dataset_id"], row["category"], row["local_path"], dest.stat().st_size,
                    row["md5"], time.strftime("%Y-%m-%dT%H:%M:%S")])


def load_plan(category="", only=""):
    with open(PLAN, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if category:
        rows = [r for r in rows if r["category"] == category]
    if only:
        want = {k.strip() for k in only.split(",") if k.strip()}
        rows = [r for r in rows if r["dataset_id"] in want]
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--category", default="")
    ap.add_argument("--only", default="")
    ap.add_argument("--connections", type=int, default=8)
    ap.add_argument("--chunk-mb", type=int, default=2)
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()

    rows = load_plan(a.category, a.only)
    if a.list:
        for r in rows:
            mb = int(r["bytes"]) / 1e6 if r["bytes"] else 0
            present = "present" if (ROOT / r["local_path"]).exists() else "missing"
            print(f"[{r['category']:5s}] {r['dataset_id']:12s} {mb:8.1f} MB {present:8s} {r['local_path']}")
            print(f"         {r['url']}")
        return 0
    if not rows:
        log("no plan rows selected")
        return 2
    ok = [r["dataset_id"] for r in rows if download(r, a.connections, a.chunk_mb)]
    bad = [r["dataset_id"] for r in rows if r["dataset_id"] not in ok]
    log(f"\ndownloaded {len(ok)}/{len(rows)}")
    if bad:
        log("incomplete: " + ", ".join(bad))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())

