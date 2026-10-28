"""Verify files listed in results/tables/download_plan.csv that were obtained by hand
(browser, download manager, another machine) and record them.

For every plan row whose destination file exists: compare the exact byte size, the md5
when the source publishes one, and read the gzip stream to the end (catches truncated
downloads and CRC errors). Present files are logged to data/raw/download_manifest.csv;
absent files are listed as still needed. Files are never deleted or overwritten here.

Usage:
    python analysis/scripts/verify_downloads.py --list
    python analysis/scripts/verify_downloads.py
    python analysis/scripts/verify_downloads.py --only gwas
"""
import argparse
import csv
import gzip
import hashlib
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PLAN = ROOT / "results" / "tables" / "download_plan.csv"
MANIFEST = ROOT / "data" / "raw" / "download_manifest.csv"


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
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 22):
                pass
        return True, ""
    except Exception as e:
        return False, f"{type(e).__name__}"
    except OSError as e:
        return False, f"OSError: {e}"


def load_plan():
    with open(PLAN, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def check(row):
    if (row.get("note") or "").startswith("UNAVAILABLE"):
        # the repository does not publish this file at all; not a download failure
        return "n/a", "not published upstream (see note in the plan)"
    dest = ROOT / row["local_path"]
    need = int(row["bytes"]) if row.get("bytes") else None
    md5 = (row.get("md5") or "").strip()
    if not dest.exists():
        return "missing", f"{row['local_path']} ({need / 1e6:.1f} MB)" if need else row["local_path"]
    got = dest.stat().st_size
    if need and got != need:
        return "bad-size", f"{got} != {need} bytes ({got - need:+d})"
    if md5 and md5 != "-":
        g = md5_of(dest)
        if g != md5:
            return "bad-md5", f"{g} != {md5}"
    if dest.suffix in (".gz", ".tgz"):
        ok, msg = gzip_ok(dest)
        if not ok:
            return "bad-gzip", msg
    return "ok", md5_of(dest) if not md5 or md5 == "-" else md5


def append_manifest(rows):
    new = not MANIFEST.exists()
    with open(MANIFEST, "a", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["dataset_id", "category", "local_path", "bytes", "md5", "verified_at"])
        for r in rows:
            dest = ROOT / r["local_path"]
            w.writerow([r["dataset_id"], r["category"], r["local_path"], dest.stat().st_size,
                        md5_of(dest), time.strftime("%Y-%m-%dT%H:%M:%S")])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--only", default="", help="filter by category, e.g. gwas or geo")
    a = ap.parse_args()

    if not PLAN.exists():
        print(f"plan not found: {PLAN.relative_to(ROOT).as_posix()}")
        return 2
    plan = [r for r in load_plan() if not a.only or r["category"] == a.only]

    if a.list:
        for r in plan:
            mb = int(r["bytes"]) / 1e6 if r.get("bytes") else 0
            print(f"[{r['category']:5s}] {r['dataset_id']:12s} {mb:8.1f} MB  -> {r['local_path']}")
            print(f"         {r['url']}")
            if r.get("md5") and r["md5"] != "-":
                print(f"         md5 {r['md5']}")
        tot = sum(int(r["bytes"]) for r in plan if r.get("bytes"))
        print(f"\n{len(plan)} files, {tot / 1e9:.2f} GB total")
        return 0

    good, bad, missing, na = [], [], [], []
    for r in plan:
        status, detail = check(r)
        line = f"{r['dataset_id']:12s} {status.upper():9s} {detail}"
        print(line)
        if status == "ok":
            good.append(r)
        elif status == "n/a":
            na.append(r)
        elif status == "missing":
            missing.append(r)
        else:
            bad.append((r, detail))

    if good:
        append_manifest(good)
    print(f"\nverified {len(good)}/{len(plan)}; corrupt {len(bad)}; missing {len(missing)}; "
          f"not published {len(na)}")
    if bad:
        print("corrupt (delete and re-download): " + ", ".join(r["dataset_id"] for r, _ in bad))
    if missing:
        print("still needed: " + ", ".join(r["dataset_id"] for r in missing))
        print(f"see the URLs with: python analysis/scripts/verify_downloads.py --list")
    return 0 if not (bad or missing) else 1


if __name__ == "__main__":
    sys.exit(main())

