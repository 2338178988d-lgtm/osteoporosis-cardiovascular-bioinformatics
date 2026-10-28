"""Module 10 鈥?why is the FinnGen effect twice the HERMES effect?

The eBMD -> HF estimate is beta = -0.034 in HERMES (Levin 2022, UKB-heavy) and
-0.064 in FinnGen R12 (no UKB). Two competing explanations:

  (i)  systematic scaling 鈥?the whole instrument-effect distribution is shifted,
       as expected from a difference in case definition or population;
  (ii) locus-specific heterogeneity 鈥?a handful of variants behave differently.

If (i), the per-instrument effects should be tightly correlated with a slope
well above 1 and no leverage points.  If (ii), the slope should be unstable and
a few instruments should carry the difference.

Outputs (results/tables/10-concordance/):
  instrument_effect_concordance.csv   per-instrument effects in both cohorts
  concordance_summary.csv             slope, correlation, leverage diagnostics
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "results" / "tables"
OUT = TAB / "10-concordance"


def paired_effects():
    """Per-instrument exposure/outcome effects for both heart-failure cohorts."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "mr", ROOT / "analysis" / "02-mr" / "scripts" / "03_mr.py")
    mr = importlib.util.module_from_spec(spec)
    sys.modules["mr"] = mr
    # 03_mr.py imports its sibling _sumstats, so that directory must be importable
    mr_dir = str(ROOT / "analysis" / "02-mr" / "scripts")
    if mr_dir not in sys.path:
        sys.path.insert(0, mr_dir)
    spec.loader.exec_module(mr)

    inst = mr.load_instruments("eBMD")
    hf = mr.outcome_at_snps("HF", set(inst["SNP"]))
    fg = mr.outcome_at_snps("FG_HF", set(inst["SNP"]))
    if hf is None or fg is None:
        return None

    # load_instruments already renames beta->bx, se->sx
    m = (inst[["SNP", "bx", "sx"]]
         .merge(hf.rename(columns={"beta": "b1", "se": "s1", "p": "p1"}), on="SNP")
         .merge(fg.rename(columns={"beta": "b2", "se": "s2", "p": "p2"}), on="SNP"))
    return m.dropna(subset=["bx", "sx", "b1", "s1", "b2", "s2"])


def ratio_estimates(m):
    """IVW estimate of the HF effect: Wald ratios b_y/b_x, weighted by bx^2/se_y^2."""
    out = {}
    for tag, b, s in [("HERMES", "b1", "s1"), ("FinnGen", "b2", "s2")]:
        wald = m[b] / m["bx"]
        w = m["bx"] ** 2 / m[s] ** 2
        b_ivw = float((w * wald).sum() / w.sum())
        se = float(np.sqrt(1.0 / w.sum()))
        out[tag] = (b_ivw, se)
    return out


def leverage(m):
    """Drop each instrument in turn and record how far the slope moves."""
    base = ratio_estimates(m)["FinnGen"][0]
    shifts = []
    for i in range(len(m)):
        sub = m.drop(m.index[i])
        b = ratio_estimates(sub)["FinnGen"][0]
        shifts.append(b - base)
    shifts = np.array(shifts)
    j = int(np.argmax(np.abs(shifts)))
    return base, float(shifts[j]), m.iloc[j]["SNP"], len(m)


def concordance(m):
    """Do the two cohorts rank the instrument Wald ratios the same way?

    Correlation is computed on the Wald ratio (b_y/b_x) of each instrument, which
    is the quantity that should agree if the effect is the same in both cohorts.
    The regression slope tells us how much larger FinnGen's effect is.
    """
    x = (m["b1"] / m["bx"]).to_numpy()
    y = (m["b2"] / m["bx"]).to_numpy()
    w = 1.0 / ((m["s1"] / m["bx"]) ** 2)
    slope, intercept = np.polyfit(x, y, 1, w=np.sqrt(w))
    r = float(stats.pearsonr(x, y)[0])
    rho = float(stats.spearmanr(x, y)[0])
    m = m.assign(wald_hermes=x, wald_finngen=y)
    return float(slope), float(intercept), r, rho, m


def error_corrected(m):
    """Is the modest per-instrument agreement just measurement error?

    The Wald ratio of an instrument with a tiny exposure effect has a very noisy
    HERMES estimate.  Pearson attenuation due to measurement error in x is
    r_true = r_obs / sqrt(reliability).  We estimate reliability from the two
    cohorts' own precision and also stratify by exposure-effect strength.
    """
    x = m["b1"] / m["bx"]
    y = m["b2"] / m["bx"]
    sx = m["s1"] / m["bx"].abs()
    sy = m["s2"] / m["bx"].abs()

    r_obs = float(stats.pearsonr(x, y)[0])
    var_x = float(x.var(ddof=1))
    var_y = float(y.var(ddof=1))
    # measurement-error variance in each Wald ratio
    mean_sx2 = float((sx ** 2).mean())
    mean_sy2 = float((sy ** 2).mean())
    rel_x = max(1e-9, 1.0 - mean_sx2 / var_x)
    rel_y = max(1e-9, 1.0 - mean_sy2 / var_y)
    r_true = r_obs / np.sqrt(rel_x * rel_y)

    # strongest instruments only: precision in the exposure effect
    top = m.nlargest(len(m) // 5, "bx")
    xt, yt = top["b1"] / top["bx"], top["b2"] / top["bx"]
    r_top = float(stats.pearsonr(xt, yt)[0]) if len(top) > 5 else np.nan

    return r_obs, r_true, rel_x, rel_y, r_top, len(top)


def between_cohort_test(est):
    """Is the FinnGen estimate actually different from the HERMES estimate?

    The ratio 1.77 looks large, but both estimates carry SEs.  A two-sample z
    test on the difference is the honest test of "the cohorts disagree".
    """
    b1, s1 = est["HERMES"]
    b2, s2 = est["FinnGen"]
    diff = b2 - b1
    se = float(np.hypot(s1, s2))
    z = diff / se
    p = float(2 * stats.norm.sf(abs(z)))
    return diff, se, z, p


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    m = paired_effects()
    if m is None or len(m) < 20:
        print("not enough shared instruments")
        return 1
    print(f"{len(m)} instruments shared by both heart-failure cohorts")

    m.to_csv(OUT / "instrument_effect_concordance.csv", index=False,
             float_format="%.6g")

    est = ratio_estimates(m)
    print(f"\n  HERMES  IVW b={est['HERMES'][0]:+.4f} (se {est['HERMES'][1]:.4f})")
    print(f"  FinnGen IVW b={est['FinnGen'][0]:+.4f} (se {est['FinnGen'][1]:.4f})")
    ratio = est["FinnGen"][0] / est["HERMES"][0]
    print(f"  ratio FinnGen/HERMES = {ratio:.2f}")

    diff, se_d, z, p_diff = between_cohort_test(est)
    print(f"\n  difference = {diff:+.4f} (se {se_d:.4f}), z = {z:+.2f}, p = {p_diff:.3f}")

    slope, intercept, r, rho, m = concordance(m)
    print(f"\n  precision-weighted slope = {slope:.3f} (intercept {intercept:+.4f})")
    print(f"  Pearson r = {r:.3f}   Spearman rho = {rho:.3f}")

    r_obs, r_true, rel_x, rel_y, r_top, n_top = error_corrected(m)
    print(f"\n  agreement corrected for measurement error:")
    print(f"    observed r = {r_obs:.3f}  -> disattenuated r = {r_true:.3f}")
    print(f"    reliability: HERMES {rel_x:.3f}, FinnGen {rel_y:.3f}")
    print(f"    top {n_top} strongest instruments: r = {r_top:.3f}")

    base, worst, worst_snp, n = leverage(m)
    print(f"\n  largest single-instrument shift: {worst:+.4f} from {worst_snp}")
    print(f"  effect without it: {base + worst:+.4f}  (relative "
          f"{100 * abs(worst / base):.1f} %)")

    row = {
        "n_instruments": len(m),
        "hermes_ivw": est["HERMES"][0],
        "finngen_ivw": est["FinnGen"][0],
        "ratio": ratio,
        "between_cohort_diff": diff,
        "between_cohort_se": se_d,
        "between_cohort_z": z,
        "between_cohort_p": p_diff,
        "slope_finngen_on_hermes": slope,
        "intercept": intercept,
        "pearson_r": r,
        "spearman_rho": rho,
        "pearson_r_disattenuated": r_true,
        "reliability_hermes": rel_x,
        "reliability_finngen": rel_y,
        "pearson_r_top_quintile": r_top,
        "n_top_quintile": n_top,
        "max_single_loco_shift": worst,
        "max_shift_snp": worst_snp,
        "max_shift_relative_pct": 100 * abs(worst / base),
    }
    pd.DataFrame([row]).to_csv(OUT / "concordance_summary.csv", index=False,
                               float_format="%.6g")
    print(f"\nwrote {OUT / 'concordance_summary.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

