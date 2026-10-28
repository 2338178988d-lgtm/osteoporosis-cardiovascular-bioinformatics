"""Pooling and reporting helpers for module 15 (no external dependency).

Inverse-variance pooling with DerSimonian-Laird heterogeneity, which is all that six
per-file estimates need: the question is whether the files inside a modality agree with each
other and whether the two modalities disagree, not how to model study weights.
"""
import numpy as np
import pandas as pd
from scipy import stats


def pool(b, se):
    """Fixed- and random-effects inverse-variance pooling with Cochran's Q and I^2."""
    b = np.asarray(b, float)
    se = np.asarray(se, float)
    w = 1 / se ** 2
    b_fixed = float(np.sum(w * b) / np.sum(w))
    se_fixed = float(np.sqrt(1 / np.sum(w)))
    k = len(b)
    q = float(np.sum(w * (b - b_fixed) ** 2))
    dfree = k - 1
    q_p = float(stats.chi2.sf(q, dfree)) if dfree > 0 else np.nan
    i2 = max(0.0, (q - dfree) / q) if q > 0 and dfree > 0 else 0.0
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (q - dfree) / c) if c > 0 and dfree > 0 else 0.0
    w_star = 1 / (se ** 2 + tau2)
    b_re = float(np.sum(w_star * b) / np.sum(w_star))
    se_re = float(np.sqrt(1 / np.sum(w_star)))
    return dict(k=k, b=b_re, se=se_re, p=float(2 * stats.norm.sf(abs(b_re / se_re))),
                b_fixed=b_fixed, se_fixed=se_fixed, q=q, q_p=q_p, i2=float(i2),
                tau2=float(tau2),
                lo=float(b_re - 1.96 * se_re), hi=float(b_re + 1.96 * se_re))


def contrast(a, b):
    """Difference between two pooled estimates, and its p-value."""
    d = a["b"] - b["b"]
    se = float(np.sqrt(a["se"] ** 2 + b["se"] ** 2))
    return dict(diff=float(d), se=se, p=float(2 * stats.norm.sf(abs(d / se))),
                lo=float(d - 1.96 * se), hi=float(d + 1.96 * se))


def pools_frame(df, outcome="FG_HF"):
    """One row per modality pool plus the between-modality contrast, as a table."""
    rows, pools = [], {}
    for modality in ("ultrasound", "dxa"):
        d = df[(df["modality"] == modality) & (df["outcome"] == outcome)]
        if d.empty:
            continue
        p = pool(d["b"], d["se"])
        pools[modality] = p
        rows.append(dict(outcome=outcome, group=modality,
                         **{k: p[k] for k in ("k", "b", "se", "p", "lo", "hi", "i2", "q_p")}))
    if len(pools) == 2:
        c = contrast(pools["dxa"], pools["ultrasound"])
        rows.append(dict(outcome=outcome, group="dxa_minus_ultrasound", k=pools["dxa"]["k"]
                         + pools["ultrasound"]["k"], b=c["diff"], se=c["se"], p=c["p"],
                         lo=c["lo"], hi=c["hi"], i2=np.nan, q_p=np.nan))
    return pd.DataFrame(rows), pools


def report(df, outcome="FG_HF"):
    """Print the within-modality pools, the between-modality contrast, and the verdict."""
    print(f"\n=== pooling within modality ({outcome})")
    table, pools = pools_frame(df, outcome)
    for _, r in table.iterrows():
        print(f"  {r['group']:22s} k={int(r['k'])} beta={r['b']:+.4f} "
              f"({r['lo']:+.4f}, {r['hi']:+.4f}) p={r['p']:.3g}"
              + ("" if np.isnan(r["i2"]) else f" I2={r['i2']:.2f} Q p={r['q_p']:.3g}"))
    d = df[df["outcome"] == outcome]
    for _, r in d.iterrows():
        print(f"      {r['exposure']:8s} {r['source']:14s} b={r['b']:+.4f} "
              f"p={r['p']:.3g} power={r['power_ref']:.0%}")

    if len(pools) == 2:
        c = contrast(pools["dxa"], pools["ultrasound"])
        print("\n--- decision rules")
        us, dx = pools["ultrasound"], pools["dxa"]
        us_all_neg = bool((df[(df["modality"] == "ultrasound")
                              & (df["outcome"] == outcome)]["b"] < 0).all())
        dx_all_pos = bool((df[(df["modality"] == "dxa")
                              & (df["outcome"] == outcome)]["b"] > 0).all())
        if us_all_neg and dx_all_pos and c["p"] < 0.05 and us["i2"] < 0.5 and dx["i2"] < 0.5:
            print("  RULE 1 - every ultrasound file protective, every DXA file positive, "
                  "within-modality heterogeneity low, contrast significant:")
            print("           the effect belongs to the phenotype, not to one file. "
                  "'Heel-ultrasound bone density' survives and the contrast is a finding.")
        elif not us_all_neg and not dx_all_pos:
            print("  RULE 2 - at least one modality disagrees with itself: source dominates "
                  "modality, and the exposure-specific wording must be withdrawn.")
        elif dx["lo"] <= 0 <= dx["hi"]:
            print("  RULE 3 - the DXA pool includes zero: the reversal is not established; "
                  "report 'DXA did not show protection' and demote the modality claim.")
        else:
            print("  no rule fires cleanly - read the table and state which one is closest")
    return table


def report_concordance(df):
    """Reciprocal prediction: files of the same phenotype must predict each other."""
    print("\n=== reciprocal prediction between exposures")
    for _, r in df.iterrows():
        flag = "" if r["b"] > 0.1 else "   <-- weak: check this file"
        print(f"  {r['pair_label']:20s} b={r['b']:+.3f} p={r['p']:.3g} "
              f"nsnp={r['nsnp']}{flag}")

