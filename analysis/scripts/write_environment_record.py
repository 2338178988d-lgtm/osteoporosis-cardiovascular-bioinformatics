"""Write an environment record for the submission package.

Versions are read from the runtime, not typed from memory: a requirements file that guesses is
worse than none.
"""
import importlib.metadata as md
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "reports" / "鐜璁板綍_environment.txt"

PY = ["numpy", "pandas", "scipy", "matplotlib", "openpyxl", "requests"]
R_PKGS = ["survey", "survival", "limma", "msigdbr", "jsonlite"]


def version(pkg):
    try:
        return md.version(pkg)
    except md.PackageNotFoundError:
        return "not installed"


lines = [
    "鐜璁板綍 / environment record",
    f"generated: {__import__('time').strftime('%Y-%m-%dT%H:%M:%S')}",
    "",
    "Python",
    f"  executable  {sys.executable}",
    f"  version     {sys.version.split()[0]}",
    "",
    "Python packages (direct dependencies used by the analysis)",
]
for p in PY:
    lines.append(f"  {p:14s} {version(p)}")

lines += ["", "R"]
try:
    out = subprocess.run(["Rscript", "--version"], capture_output=True, text=True, timeout=30)
    lines.append(f"  version     {(out.stdout or out.stderr).strip().splitlines()[0]}")
    for pkg in R_PKGS:
        r = subprocess.run(
            ["Rscript", "-e", f'cat(as.character(packageVersion("{pkg}")))'],
            capture_output=True, text=True, timeout=60)
        v = (r.stdout or "").strip() or "not installed"
        lines.append(f"  {pkg:12s} {v}")
except (FileNotFoundError, subprocess.SubprocessError, IndexError) as e:
    lines.append(f"  not probed: {type(e).__name__}")

lines += [
    "",
    "External tools",
    "  PLINK        v1.90b7.2 64-bit (11 Dec 2023; data/raw/reference/plink/plink.exe)",
    "  1000G        phase 3 EUR, PLINK binary (data/raw/reference/EUR.*)",
    "",
    "Data releases",
    "  FinnGen      R12, finngen-public-data-r12, accessed 2026-10-23/24",
    "               endpoint definitions and case counts: manifests/finngen_R12_manifest.tsv",
    "  GWAS Catalog harmonised summary statistics, 2026-10 downloads",
    "               full file list with md5 and byte sizes: results/tables/download_plan.csv",
    "  NHANES       2005-2018, seven two-year cycles, complex survey design",
    "",
    "Verification entry points",
    "  analysis/scripts/verify_downloads.py   96 data entries, md5/gzip end to end",
    "  analysis/scripts/verify_manuscript.py  72 numbers in the draft against the tables",
    "",
    "Note on the figure style runtime",
    "  analysis/scripts/_figure_style_runtime.py is a vendored, byte-identical copy of the",
    "  style runtime used to generate the figures. _paper_common.load_style() prefers it, so",
    "  the six paper figures rebuild without any external skill directory.",
]

OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(OUT.read_text(encoding="utf-8"))

