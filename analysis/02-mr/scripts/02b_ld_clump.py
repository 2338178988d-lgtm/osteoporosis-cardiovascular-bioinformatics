"""LD clumping with a local 1000G reference panel (PLINK).

02_extract_instruments.py could only prune by distance because Ensembl's LD endpoints
return one variant at a time. This script does the real thing: write the significant
variants to PLINK format, run `plink --clump`, and rewrite the instrument tables so that
each locus contributes one variant.

Input : results/tables/02-mr/instruments_<KEY>_all.csv  (written by 02 with --keep-all)
Output: results/tables/02-mr/instruments_<KEY>.csv      (LD-clumped, one row per locus)
        results/tables/02-mr/ld_clump_summary.csv

Usage:
    python analysis/02-mr/scripts/02b_ld_clump.py --datasets eBMD,CAD,MI,STROKE,HF,OP
    python analysis/02-mr/scripts/02b_ld_clump.py --datasets eBMD --r2 0.001 --kb 10000
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "results" / "tables" / "02-mr"
REF = ROOT / "data" / "raw" / "reference"
PLINK = REF / "plink" / "plink.exe"
EUR = REF / "EUR"


def write_plink_inputs(d, path):
    """Assoc file for --clump: SNP, P; the .assoc format PLINK expects.

    Chromosomes are written as numbers, and variants outside chr1-22 are omitted:
    the EUR panel has no X, so PLINK cannot clump the X-chromosome hits that some
    traits (estradiol) contribute.
    """
    on_panel = d[pd.to_numeric(d["chr"], errors="coerce").between(1, 22).to_numpy()]
    if len(on_panel) < len(d):
        print(f"  {len(d) - len(on_panel)} variants outside chr1-22 skipped")
    with open(path, "w", encoding="utf-8") as f:
        f.write("CHR SNP BP P\n")
        for _, r in on_panel.iterrows():
            f.write(f"{int(float(r['chr']))} {r['SNP']} {int(r['pos'])} {r['p']:.6g}\n")


def run_clump(assoc, out_prefix, r2, kb, threads=4):
    cmd = [str(PLINK), "--bfile", str(EUR), "--clump", str(assoc),
           "--clump-p1", "1", "--clump-r2", str(r2), "--clump-kb", str(kb),
           "--clump-snp-field", "SNP", "--clump-field", "P",
           "--out", str(out_prefix), "--allow-no-sex", "--memory", "4000"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", default="eBMD,CAD,MI,STROKE,HF,OP")
    ap.add_argument("--r2", type=float, default=0.001)
    ap.add_argument("--kb", type=int, default=10000)
    a = ap.parse_args()

    if not PLINK.exists():
        print(f"plink not found at {PLINK}")
        return 2
    if not Path(f"{EUR}.bed").exists():
        print(f"reference panel not found at {EUR}")
        return 2

    summary = []
    for key in [k.strip() for k in a.datasets.split(",") if k.strip()]:
        src = TAB / f"instruments_{key}_all.csv"
        if not src.exists():
            print(f"{key}: {src.name} missing (run 02_extract_instruments.py --keep-all)")
            continue
        d = pd.read_csv(src)
        print(f"\n=== {key}: {len(d)} candidates")
        with tempfile.TemporaryDirectory() as tmp:
            assoc = Path(tmp) / f"{key}.assoc"
            write_plink_inputs(d, assoc)
            out = Path(tmp) / f"{key}_clumped"
            res = run_clump(assoc, out, a.r2, a.kb)
            clumped = Path(f"{out}.clumped")
            if not clumped.exists():
                print("  plink failed: " + (res.stdout or "")[-1500:] + (res.stderr or "")[-800:])
                continue
            c = pd.read_csv(clumped, sep=r"\s+", engine="python")
            if "SNP" not in c.columns:
                print(f"  unexpected plink output columns: {list(c.columns)}")
                continue
            keep = set(c["SNP"].astype(str))
            inst = d[d["SNP"].isin(keep)].sort_values("p").reset_index(drop=True)
            n_before = len(d)
            # the analysis scripts read instruments_<KEY>.csv: replace it, keep the
            # unclumped candidate list alongside for reproducibility
            inst.to_csv(TAB / f"instruments_{key}.csv", index=False)
        summary.append(dict(dataset=key, candidates=n_before, clumped=len(inst),
                            removed=n_before - len(inst), r2=a.r2, kb=a.kb,
                            min_F=round(inst["F"].min(), 1),
                            sum_r2=round(inst["r2"].sum(), 4)))
        print(f"  r2<{a.r2} within {a.kb} kb -> {len(inst)} independent instruments "
              f"(removed {n_before - len(inst)}), min F = {inst['F'].min():.1f}, "
              f"Sum R2 = {inst['r2'].sum():.3f}")
        if (TAB / "instruments_summary.csv").exists():
            pass
    if summary:
        s = pd.DataFrame(summary)
        s.to_csv(TAB / "ld_clump_summary.csv", index=False)
        print("\n" + s.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())

