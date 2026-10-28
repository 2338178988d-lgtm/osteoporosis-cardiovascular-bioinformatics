"""Shared estimator wrapper for module 13: one pair, IVW + Egger + empirical power.

Every estimate carries the power it had to detect the reference effect of |b| = 0.064
(the FinnGen eBMD estimate) and the smallest effect it could detect at 80 % power, because
a null without those two numbers cannot be read.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MR_DIR = HERE.parents[1] / "02-mr" / "scripts"
sys.path.insert(0, str(MR_DIR))
_spec = importlib.util.spec_from_file_location("mr", MR_DIR / "03_mr.py")
mr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mr)

import _power  # noqa: E402


def have(key):
    import _sumstats
    return key in _sumstats.FILES and _sumstats.path_of(key).exists()


def implied_sd(key, nmax=300_000):
    """Estimate the outcome trait's SD from its own summary statistics.

    For a quantitative trait, se^2 ~ 1/(2 f (1-f) N) when the trait is standardised; the
    ratio se / sqrt(1/(2f(1-f)N)) therefore estimates the SD in the trait's own units. This
    is what lets a raw beta (grams of LV mass) be reported in SD units.
    """
    import gzip
    import _sumstats
    path = _sumstats.path_of(key)
    fallback = _sumstats.TRAIT_N[key][0]
    se_l, f_l, n_l = [], [], []
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        c_se, c_f = idx["standard_error"], idx["effect_allele_frequency"]
        c_n = idx.get("n") or idx.get("N")
        for k, line in enumerate(f):
            if k >= nmax:
                break
            p = line.rstrip("\n").split("\t")
            try:
                se, fr = float(p[c_se]), float(p[c_f])
                n = float(p[c_n]) if c_n else float(fallback)
            except (ValueError, IndexError, TypeError):
                continue
            if se <= 0 or not (0.01 < fr < 0.99) or n <= 0:
                continue
            se_l.append(se)
            f_l.append(fr)
            n_l.append(n)
    if not se_l:
        return np.nan
    se = np.array(se_l)
    fr = np.array(f_l)
    nn = np.array(n_l)
    return float(np.median(se * np.sqrt(2 * fr * (1 - fr) * nn)))


def estimate(exp_key, out_key, inst=None, label="", min_snps=8, reference=0.064, out=None):
    """IVW (random effects) + MR-Egger for one exposure/outcome pair.

    `out` overrides the local-file lookup, for outcomes that are read through the remote
    tabix reader (module 14) instead of a downloaded summary-statistics file.
    """
    if out is None and not have(out_key):
        print(f"  {out_key}: file absent, skipped")
        return None
    if inst is None:
        inst = mr.load_instruments(exp_key)
    if inst is None or inst.empty:
        print(f"  {exp_key}: no instruments")
        return None
    out = out if out is not None else mr.outcome_at_snps(out_key, inst["SNP"].tolist())
    if out is None or out.empty:
        print(f"  {out_key}: no outcome matches")
        return None
    d = mr.harmonise(inst, out)
    if len(d) < min_snps:
        print(f"  {out_key}: only {len(d)} usable SNPs")
        return None
    bx, by = d["bx"].to_numpy(float), d["by"].to_numpy(float)
    sx, sy = d["sx"].to_numpy(float), d["sy"].to_numpy(float)
    r = mr.ivw(bx, by, sx, sy, random=True)
    e = mr.egger(bx, by, sx, sy)
    se = float(r["se"])
    row = dict(pair=label or f"{exp_key}->{out_key}", exposure=exp_key, outcome=out_key,
               nsnp=len(d), b=float(r["b"]), se=se, p=float(r["p"]),
               lo=float(r["b"]) - 1.96 * se, hi=float(r["b"]) + 1.96 * se,
               q_p=float(r["q_p"]), egger_intercept_p=float(e["intercept_p"]),
               power_ref=_power.power(se, reference),
               mde_80=_power.mde(se, 0.80))
    print(f"  {row['pair']:36s} n={len(d):4d}  b={row['b']:+.4f} "
          f"({row['lo']:+.4f}, {row['hi']:+.4f})  p={row['p']:.3g}  "
          f"power={row['power_ref']:.0%}  MDE80={row['mde_80']:.3f}")
    return row

