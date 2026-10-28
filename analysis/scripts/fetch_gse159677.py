"""Stream GSE159677_RAW.tar and keep only the six per-sample filtered matrices.

The per-sample bundles listed in GEO's filelist.txt are not individually downloadable
(each returns 404); they exist only inside the 2.1 GB RAW.tar. Downloading the whole
archive and unpacking it would put ~2.1 GB of raw sequencing bundles on disk of which we
need 340 MB, and most of the archive is raw molecule-info HDF5 files we never use.

So the archive is fetched as a stream and passed straight to tarfile, which extracts the
six `*_featurebcmatrixfiltered.tar.gz` members and skips everything else. Nothing but the
six bundles is written.

Output: data/raw/scrna/<GSM>/<original name>  (one bundle per sample)
        data/raw/scrna/GSM_samples.csv         (sample -> patient / tissue)

Usage:
    python analysis/scripts/fetch_gse159677.py
    python analysis/scripts/fetch_gse159677.py --list
"""
import argparse
import csv
import sys
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "data" / "raw" / "scrna"
URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159677/suppl/"
       "GSE159677_RAW.tar")
UA = {"User-Agent": "Mozilla/5.0 (research; module 07-scrna)"}

# from GSM4837523..GSM4837528; AC = atherosclerotic (calcified) core, PA = proximal adjacent
SAMPLES = {
    "GSM4837523": ("Patient 1", "AC"),
    "GSM4837524": ("Patient 1", "PA"),
    "GSM4837525": ("Patient 2", "AC"),
    "GSM4837526": ("Patient 2", "PA"),
    "GSM4837527": ("Patient 3", "AC"),
    "GSM4837528": ("Patient 3", "PA"),
}


class HttpFile:
    """Minimal seekable-free file object over an HTTP response, for tarfile streaming."""

    def __init__(self, url, chunk=1 << 20):
        self._resp = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300)
        self._chunk = chunk
        self._buf = b""
        self.bytes_read = 0
        self.total = int(self._resp.headers.get("Content-Length") or 0)

    def read(self, n=-1):
        while n < 0 or len(self._buf) < n:
            piece = self._resp.read(self._chunk if n < 0 else max(n - len(self._buf), 1 << 20))
            if not piece:
                break
            self._buf += piece
            self.bytes_read += len(piece)
        if n < 0:
            out, self._buf = self._buf, b""
            return out
        out, self._buf = self._buf[:n], self._buf[n:]
        return out

    def close(self):
        self._resp.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()

    if a.list:
        for gsm, (pat, tis) in sorted(SAMPLES.items()):
            have = list((DEST / gsm).glob("*featurebcmatrixfiltered.tar.gz")) if (DEST / gsm).exists() else []
            print(f"{gsm}  {pat:9s} {tis}  {'present' if have else 'missing'}")
        return 0

    DEST.mkdir(parents=True, exist_ok=True)
    keep_total, seen = 0, []
    stream = HttpFile(URL)
    print(f"streaming {stream.total / 1e9:.2f} GB; keeping only the six filtered matrices")
    try:
        with tarfile.open(fileobj=stream, mode="r|") as tf:
            for member in tf:
                name = Path(member.name).name
                if not name.endswith("_featurebcmatrixfiltered.tar.gz"):
                    continue
                gsm = name[:10]
                if gsm not in SAMPLES:
                    print(f"  skipping unrecognised member {name}")
                    continue
                outdir = DEST / gsm
                outdir.mkdir(parents=True, exist_ok=True)
                target = outdir / name
                fh = tf.extractfile(member)
                if fh is None:
                    continue
                with open(target, "wb") as out:
                    while True:
                        block = fh.read(1 << 20)
                        if not block:
                            break
                        out.write(block)
                keep_total += target.stat().st_size
                seen.append(gsm)
                print(f"  {gsm} {name[:44]}... {target.stat().st_size / 1e6:.1f} MB "
                      f"(streamed {stream.bytes_read / 1e9:.2f} GB)")
                if len(seen) == len(SAMPLES):
                    print("  all six found; closing the stream early")
                    break
    finally:
        stream.close()

    with open(DEST / "GSM_samples.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gsm", "patient", "tissue"])
        for gsm in sorted(SAMPLES):
            w.writerow([gsm, SAMPLES[gsm][0], SAMPLES[gsm][1]])

    missing = [g for g in SAMPLES if g not in seen]
    print(f"\nkept {len(seen)}/{len(SAMPLES)} bundles, {keep_total / 1e6:.1f} MB on disk")
    if missing:
        print(f"missing: {', '.join(missing)}")
    return 0 if not missing else 1


if __name__ == "__main__":
    sys.exit(main())

