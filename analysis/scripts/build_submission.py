"""Build the submission package.

A submission package is not a backup: it contains what a reader or reviewer needs to check the
work, and omits what they can refetch. So:

  included   manuscript, figures, every result table carrying a number in the paper, the code
             that produced them (with READMEs and self-tests), the reports, the data manifest
             with per-file checksums, and the verification scripts
  excluded   data/raw (~40 GB of public GWAS/GEO files, re-downloadable through the manifest),
             data/processed (derived caches), .venv, .cache, .wisp, runs
  excluded   results/tables/*_all.csv - the unclumped candidate lists. They are intermediates,
             not results: every instrument set used in the paper is the clumped
             instruments_<key>.csv

Layout of the package:

  README.md                     the package summary
  MANIFEST_SHA256.txt           per-file checksums
  manuscript/                   the draft and its writing notes
  figures/paper/                the six paper figures
  figures/                      the twelve process figures they were distilled from
  results/tables/, reports/     every table and report
  analysis/                     the code

Usage:
  python analysis/scripts/build_submission.py
  python analysis/scripts/build_submission.py --keep-staging   # leave submission/ for review
"""
import argparse
import hashlib
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STAGE = ROOT / "submission"
OUT = ROOT / "results" / "submission"
NAME = "骨密度心衰研究_投稿包_v1"

EXCLUDE_PARTS = {".venv", ".cache", ".wisp", "runs", "__pycache__", ".git",
                 "submission", ".downloads", "_tbi", ".chunks", "_replaced", "backup"}
EXCLUDE_SUFFIX = {".part", ".pyc", ".assembling", ".corrupt", ".tbi"}


def keep(rel: Path) -> bool:
    if any(part in EXCLUDE_PARTS for part in rel.parts):
        return False
    if rel.suffix in EXCLUDE_SUFFIX or rel.name.endswith(".log"):
        return False
    # only the unclumped candidate lists are intermediates. Matching on "_all.csv" alone also
    # caught big_cohorts_all.csv, which is a result table the manuscript and fig 2 both use.
    if rel.name.startswith("instruments_") and rel.name.endswith("_all.csv"):
        return False
    return True


def collect() -> list[tuple[Path, Path]]:
    """(source, destination-relative-to-package) pairs."""
    pairs: list[tuple[Path, Path]] = []

    def add_tree(src_root: Path, dest_root: Path):
        if not src_root.exists():
            return
        for f in src_root.rglob("*"):
            if not f.is_file():
                continue
            rel = f.relative_to(ROOT)
            if not keep(rel):
                continue
            pairs.append((f, dest_root / f.relative_to(src_root)))

    add_tree(ROOT / "analysis", Path("analysis"))
    add_tree(ROOT / "results" / "tables", Path("results/tables"))
    add_tree(ROOT / "literature", Path("literature"))

    # figures: paper figures in their own folder, process figures at the top of figures/
    for f in sorted((ROOT / "figures" / "paper").glob("*")):
        if f.is_file():
            pairs.append((f, Path("figures/paper") / f.name))
    for f in sorted((ROOT / "figures").glob("fig*")):
        if f.is_file():
            pairs.append((f, Path("figures") / f.name))

    # reports: the manuscript goes to manuscript/, the rest to results/reports/
    reports = ROOT / "results" / "reports"
    for f in sorted(reports.glob("*.md")):
        if f.name.startswith(("论文草稿_", "论文草稿.")):
            pairs.append((f, Path("manuscript") / f.name))
        elif f.name.startswith("投稿包说明"):
            continue
        else:
            pairs.append((f, Path("results/reports") / f.name))

    # the package README
    readme = reports / "投稿包说明_README.md"
    if readme.exists():
        pairs.append((readme, Path("README.md")))

    # single files worth carrying
    for extra in ["results/tables/download_plan.csv",
                  "data/processed/finngen_R12_manifest.tsv"]:
        p = ROOT / extra
        if p.exists():
            pairs.append((p, Path(extra.replace("data/processed/", "manifests/"))))

    return sorted(set(pairs), key=lambda p: str(p[1]))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep-staging", action="store_true")
    a = ap.parse_args()

    if STAGE.exists():
        shutil.rmtree(STAGE)
    STAGE.mkdir(parents=True)

    pairs = collect()
    total = 0
    for src, dest_rel in pairs:
        dest = STAGE / dest_rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        total += src.stat().st_size
    print(f"staged {len(pairs)} files, {total / 1e6:.1f} MB")

    man = STAGE / "MANIFEST_SHA256.txt"
    lines = [f"{sha256(f)}  {f.relative_to(STAGE).as_posix()}"
             for f in sorted(STAGE.rglob("*"))
             if f.is_file() and f.name != "MANIFEST_SHA256.txt"]
    man.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"manifest: {len(lines)} entries")

    OUT.mkdir(parents=True, exist_ok=True)
    zip_path = OUT / f"{NAME}.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for f in sorted(STAGE.rglob("*")):
            if f.is_file():
                zf.write(f, f.relative_to(STAGE).as_posix())
    digest = sha256(zip_path)
    (OUT / f"{NAME}.zip.sha256").write_text(f"{digest}  {NAME}.zip\n", encoding="utf-8")

    print(f"\narchive : {zip_path.relative_to(ROOT).as_posix()}")
    print(f"  size    {zip_path.stat().st_size / 1e6:.1f} MB")
    print(f"  files   {len(lines)}")
    print(f"  sha256  {digest}")

    if not a.keep_staging:
        shutil.rmtree(STAGE)
        print("  staging removed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
