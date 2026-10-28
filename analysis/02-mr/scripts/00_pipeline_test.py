"""End-to-end check of the 02-mr pipeline on simulated summary statistics.

Writes synthetic GWAS-Catalog-harmonised files with a known causal effect, runs the
instrument extraction and the MR estimators on them, and checks that the inverse-
variance weighted estimate recovers the simulated effect.

Usage: python analysis/02-mr/scripts/00_pipeline_test.py
"""
import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
TRUE_BETA = 0.30
N_INS = 60


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, HERE / fname)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make_sumstats(n_snps, seed):
    rng = np.random.default_rng(seed)
    bx = rng.normal(0, 0.05, n_snps)
    bx = np.sign(bx) * np.maximum(np.abs(bx), 0.02)      # |z| > 10, as for p < 5e-8 hits
    sx = np.full(n_snps, 0.002)
    sy = np.full(n_snps, 0.02)
    by = TRUE_BETA * bx + rng.normal(0, sy)
    return bx, sx, by, sy


def write_file(path, b, se, n=200000):
    df = pd.DataFrame({
        "hm_rsid": [f"rs{100000 + i}" for i in range(len(b))],
        "hm_chrom": np.repeat(1, len(b)),
        "hm_pos": np.arange(len(b)) * 2_000_000 + 1_000_000,
        "hm_effect_allele": "A",
        "hm_other_allele": "G",
        "hm_effect_allele_frequency": 0.3,
        "hm_beta": b,
        "standard_error": se,
        "p_value": 1e-12,
        "n": n,
        "beta": b,
        "se": se,
        "p": 1e-3,
    })
    df.to_csv(path, sep="\t", index=False)


def main():
    ext = load("ext", "02_extract_instruments.py")
    mr = load("mr", "03_mr.py")
    tmp = Path(tempfile.mkdtemp(prefix="mr_selftest_"))
    bx, sx, by, sy = make_sumstats(N_INS, 3)
    write_file(tmp / "GCST_TEST_eBMD.h.tsv.gz", bx, sx)
    # outcome file: same SNPs, causal effect of TRUE_BETA
    write_file(tmp / "GCST_TEST_CAD.h.tsv.gz", by, sy)

    ext.RAW, ext.FILES, ext.OUT = tmp, {"eBMD": "GCST_TEST_eBMD.h.tsv.gz"}, tmp / "tables"
    mr.RAW, mr.TAB, mr.FILES = tmp, tmp / "tables", {"eBMD": "GCST_TEST_eBMD.h.tsv.gz",
                                                     "CAD": "GCST_TEST_CAD.h.tsv.gz"}
    mr.CACHE = tmp / "cache"
    # the column reader is shared, so point it at the synthetic files too
    ext._sumstats.RAW = tmp
    ext._sumstats.FILES = {"eBMD": "GCST_TEST_eBMD.h.tsv.gz", "CAD": "GCST_TEST_CAD.h.tsv.gz"}
    sys.argv = ["x", "--datasets", "eBMD", "--clump", "distance", "--kb", "500"]
    ext.main()

    inst = pd.read_csv(tmp / "tables" / "instruments_eBMD.csv")
    assert len(inst) == N_INS, f"expected {N_INS} instruments, got {len(inst)}"

    res, s = mr.run_pair("eBMD", "CAD", verbose=True)
    ivw = res[res["method"] == "IVW (random)"].iloc[0]
    print(f"\nsimulated beta={TRUE_BETA}: this dataset gave {ivw['b']:.3f} "
          f"(se {ivw['se']:.3f}); Steiger p={s['steiger_p']:.2g}")

    # replication over independent simulated datasets: the estimator must be unbiased
    # and its 95% interval must cover the simulated effect at roughly nominal rate
    bs, cover = [], []
    for seed in range(100):
        b2, x2, y2, s2 = make_sumstats(N_INS, 100 + seed)
        r = mr.ivw(b2, y2, x2, s2)
        bs.append(r["b"])
        cover.append(abs(r["b"] - TRUE_BETA) <= 1.96 * r["se"])
    bias = float(np.mean(bs) - TRUE_BETA)
    z = bias / (np.std(bs, ddof=1) / np.sqrt(len(bs)))
    cov = float(np.mean(cover))
    print(f"{len(bs)} simulated datasets: mean estimate {np.mean(bs):.3f} "
          f"(bias {bias:+.3f}, z={z:+.2f}), 95% CI coverage {cov:.0%}")
    ok = abs(z) < 3 and cov >= 0.85
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

