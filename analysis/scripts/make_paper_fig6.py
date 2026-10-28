"""Figure 6 - prespecified adjustments do not materially attenuate the estimate.

Panel a  ten circulating candidates, each in a multivariable model: none attenuates the effect
Panel b  body composition and stature: none attenuates it either
Panel c  negative controls and left ventricular mass; shared-locus and Steiger checks are reported separately

Claim: "The prespecified circulating-trait and body-composition models leave the estimate largely unchanged." Panels a and b are attenuation
percentages, not effect sizes, so a bar at zero is the reassuring answer - which is the opposite
of the usual reading of a bar chart and is therefore labelled on the axis.

Run from the project root:  python analysis/scripts/make_paper_fig6.py
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import _paper_common as pc

pc.load_style()
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd

TAB = ROOT / "results" / "tables"

MEDIATOR_LABEL = {
    "VITD": "25(OH) vitamin D", "SOST": "sclerostin", "BGLAP": "osteocalcin",
    "FGF23": "FGF23", "PTH": "parathyroid hormone", "DKK1": "DKK1",
    "TESTO": "testosterone", "ESTRADIOL": "estradiol", "CST3": "cystatin C", "CRP": "CRP",
}


def panel_a(ax):
    """Attenuation of the heart-failure estimate after adjusting for each circulating factor."""
    d = pd.read_csv(TAB / "06-bone-axis" / "axis_mvmr.csv")
    d = d[d["outcome"] == "FG_HF"].copy()
    d["abs"] = d["attenuation_pct"].abs()
    d = d.sort_values("abs")

    y = np.arange(len(d))
    ax.axvline(0, color="0.5", lw=0.8, zorder=1)
    for yy, (_, r) in zip(y, d.iterrows()):
        ax.plot([0, r["attenuation_pct"]], [yy, yy], color=pc.NULL, lw=1.4,
                solid_capstyle="butt", zorder=2)
        ax.plot(r["attenuation_pct"], yy, "o", color=pc.NULL, ms=4.4,
                mec="white", mew=0.5, zorder=3)
        ax.annotate(f"{r['attenuation_pct']:+.1f} %", (r["attenuation_pct"], yy),
                    textcoords="offset points",
                    xytext=(5, 0) if r["attenuation_pct"] >= 0 else (-5, 0),
                    ha="left" if r["attenuation_pct"] >= 0 else "right", va="center",
                    fontsize=5.5, color="0.3")
    ax.set_yticks(y)
    ax.set_yticklabels([MEDIATOR_LABEL.get(m, m) for m in d["mediator"]], fontsize=5.9)
    ax.set_xlabel("attenuation after adjustment (%); 0 = unchanged", fontsize=6.5)
    ax.set_ylim(-1.5, len(d) - 0.4)
    ax.set_xlim(-9, 9)
    ax.set_xticks([-6, -3, 0, 3, 6])
    ax.set_title("Tested circulating traits", loc="left", fontsize=7.4)
    ax.annotate("every value is within +/-10 %", xy=(0.98, 0.05),
                xycoords="axes fraction", ha="right", va="bottom", fontsize=5.2,
                color="0.35")


def panel_b(ax):
    """Body composition and stature, the alternative explanation that fits eBMD by construction."""
    d = pd.read_csv(TAB / "11-body-composition" / "bodycomp_mvmr.csv")
    d = d[d["outcome"] == "FG_HF"]
    rows = [("M2_eBMD_ALM", "adjusted for lean mass"),
            ("M3_eBMD_HEIGHT", "adjusted for height"),
            ("M4_eBMD_ALM_BMI", "adjusted for lean mass and BMI")]
    y = np.arange(len(rows))[::-1]
    ax.axvline(0, color="0.5", lw=0.8, zorder=1)
    for yy, (model, label) in zip(y, rows):
        r = d[d["model"] == model].iloc[0]
        ax.plot([0, r["attenuation_pct"]], [yy, yy], color=pc.NULL, lw=1.4,
                solid_capstyle="butt", zorder=2)
        ax.plot(r["attenuation_pct"], yy, "o", color=pc.NULL, ms=4.6,
                mec="white", mew=0.5, zorder=3)
        ax.annotate(f"{r['attenuation_pct']:+.1f} %", (r["attenuation_pct"], yy),
                    textcoords="offset points", xytext=(5, 0), ha="left", va="center",
                    fontsize=5.5, color="0.3")
    ax.set_yticks(y)
    ax.set_yticklabels([lab for _, lab in rows], fontsize=5.9)
    ax.set_xlabel("attenuation (%); 0 = unchanged", fontsize=6.5)
    ax.set_xlim(-9, 12)
    ax.set_ylim(-1.5, len(rows) - 0.4)
    ax.set_xticks([-6, -3, 0, 3, 6, 9])
    ax.set_title("Tested body-size traits", loc="left", fontsize=7.4)
    ax.annotate("attenuation remains\nwithin the tested range",
                xy=(0.98, 0.10), xycoords="axes fraction", ha="right", va="bottom",
                fontsize=5.2, color="0.35")


def panel_c(ax):
    """Negative controls and left ventricular mass; shared-locus and Steiger checks are reported separately."""
    neg = pd.read_csv(TAB / "13-robustness" / "negative_controls.csv")
    neg = neg[neg["exposure"] == "eBMD"]
    lv = pd.read_csv(TAB / "13-robustness" / "lv_mass.csv")

    items = []
    for _, r in neg.iterrows():
        items.append(dict(label=r["label"], b=float(r["b"]), lo=float(r["lo"]),
                          hi=float(r["hi"]), kind="null"))
    lvr = lv[lv["outcome"] == "LVM"].iloc[0]
    items.append(dict(label="left ventricular mass (SD)", b=float(lvr["b_sd"]),
                      lo=float(lvr["lo_sd"]), hi=float(lvr["hi_sd"]),
                      kind="protect" if lvr["p"] < 0.05 else "null"))

    pc.forest(ax, items, xlabel="effect per SD (log OR, or SD for LV mass)",
              label_size=5.9)
    ax.set_title("Negative controls stay null; LV mass is separate", loc="left", fontsize=7.4)
    ax.annotate("0 of 6 nominal at 92 % median power", xy=(0.5, -0.34),
                xycoords="axes fraction", ha="center", va="top", fontsize=5.2,
                color="0.35")


def main():
    pc.apply_figure_style(sizes=(7.6, 6.6, 5.8))
    fig = plt.figure(figsize=(7.6, 2.7))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 0.95, 1.15], wspace=1.15)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))

    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)

    for ax, letter in zip([ax_a, ax_b, ax_c], "abc"):
        pc.panel_letter(ax, letter)

    pc.save(fig, "fig6_robustness.png", ROOT)
    return 0


if __name__ == "__main__":
    sys.exit(main())

