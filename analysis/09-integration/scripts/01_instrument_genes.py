"""Module 09 鈥?connect the genetic hits to the cells that express them.

Until now the three layers of this project never touched: the MR layer knows *which
variants* move bone density, and the single-cell layer knows *which cell types* carry a
bone-like programme, but nothing has joined them. The 600 eBMD instruments have never
been assigned to genes, so "bone-like programme in macrophages" and "eBMD lowers
heart-failure risk" sit side by side without evidence that they describe the same biology.

FinnGen summary statistics carry a `nearest_genes` column on every row, so the assignment
needs no external annotation and no network access: scan the FinnGen heart-failure file
once, keep every row whose rsID is an eBMD instrument, and read off the gene.

Output: results/tables/09-integration/instrument_genes.csv
  one row per eBMD instrument: SNP, chr, pos, nearest gene, matched flag

Usage:
  python analysis/09-integration/scripts/01_instrument_genes.py
  python analysis/09-integration/scripts/01_instrument_genes.py --self-test
"""
import argparse
import gzip
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
FINNGEN = ROOT / "data" / "raw" / "finngen" / "finngen_R12_I9_HEARTFAIL.gz"
INSTR = ROOT / "results" / "tables" / "02-mr" / "instruments_eBMD.csv"
OUT = ROOT / "results" / "tables" / "09-integration"


def scan_nearest_genes(snps, path=FINNGEN):
    """Stream a FinnGen file once and return rsID -> (chr, pos, nearest_genes)."""
    want = set(snps)
    found = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx_rs = header.index("rsids")
        idx_gene = header.index("nearest_genes")
        # FinnGen writes the first column as "#chrom"; GWAS Catalog harmonised files as
        # "chr". Accept either so the scanner works on both layouts.
        idx_chr = header.index("#chrom") if "#chrom" in header else header.index("chrom")
        idx_pos = header.index("pos")
        n = 0
        for line in f:
            n += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= idx_gene:
                continue
            rs = parts[idx_rs]
            if rs in want and rs not in found:
                found[rs] = (parts[idx_chr], parts[idx_pos], parts[idx_gene])
                if len(found) == len(want):
                    break
            if n % 4_000_000 == 0:
                print(f"    scanned {n:,} rows, found {len(found)}/{len(want)}", flush=True)
    return found, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    inst = pd.read_csv(INSTR)
    print(f"eBMD instruments: {len(inst)}")
    found, scanned = scan_nearest_genes(inst["SNP"].tolist())
    print(f"scanned {scanned:,} FinnGen rows; matched {len(found)}/{len(inst)} instruments")

    rows = []
    for _, r in inst.iterrows():
        hit = found.get(r["SNP"])
        rows.append(dict(SNP=r["SNP"], chr=int(r["chr"]), pos=int(r["pos"]),
                         beta=float(r["beta"]), se=float(r["se"]), p=float(r["p"]),
                         F=float(r["F"]),
                         nearest_gene=hit[2] if hit else "",
                         matched=hit is not None))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "instrument_genes.csv", index=False)

    n_ok = int(d["matched"].sum())
    print(f"\nmatched {n_ok}/{len(d)} instruments ({100 * n_ok / len(d):.1f}%)")
    genes = d.loc[d.matched & (d.nearest_gene != ""), "nearest_gene"]
    genes = genes[~genes.str.contains(",")]
    vc = genes.value_counts()
    print(f"distinct nearest genes: {vc.shape[0]}")
    print("\nmost frequent nearest genes across instrument loci:")
    print(vc.head(20).to_string())

    known = {"ESR1", "WNT16", "CPED1", "LRP5", "SOST", "TNFSF11", "TNFRSF11B",
             "COL1A1", "SOX6", "RSPO3", "FKBP5", "GALNT3", "RUNX2", "SP7", "ALPL"}
    print(f"\nknown BMD genes among nearest genes: {sorted(set(vc.index) & known)}")
    return 0


def selftest():
    """The scanner must return only the requested rsIDs, with the right gene."""
    rows = [
        "#chrom\tpos\tref\talt\trsids\tnearest_genes\tpval\tmlogp\tbeta\tsebeta",
        "1\t100\tA\tG\trs111\tGENE_A\t0.5\t0.3\t0.01\t0.1",
        "1\t200\tC\tT\trs222\tGENE_B\t0.1\t1.0\t0.02\t0.1",
        "2\t300\tG\tA\trs333\tGENE_C\t0.9\t0.05\t0.00\t0.1",
    ]
    with tempfile.NamedTemporaryFile("w", suffix=".gz", delete=False, encoding="utf-8") as fh:
        path = Path(fh.name)
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write("\n".join(rows) + "\n")
    found, n = scan_nearest_genes(["rs111", "rs333", "rs999"], path=path)
    path.unlink(missing_ok=True)
    print(f"selftest: scanned {n} rows, found {sorted(found)}")
    ok = (set(found) == {"rs111", "rs333"}
          and found["rs111"][2] == "GENE_A"
          and found["rs333"][2] == "GENE_C")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

