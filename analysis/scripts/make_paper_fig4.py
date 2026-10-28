"""Figure 4 - the exposure question: is this a bone-density effect, or an ultrasound one?

Panel a  six bone-density GWAS files against FinnGen heart failure, three per modality, with the
         within-modality pool drawn as a diamond
Panel b  the same six against HERMES, where the DXA group disagrees with itself
Panel c  reciprocal prediction between the files: within a modality they predict each other at
         about 1 SD per SD, across modalities at 0.45, so the comparison in a and b is between
         two measurements of one phenotype rather than between unrelated traits

Claim: "Every ultrasound file is protective; the DXA files are neither consistently opposite nor
consistently null." Panels a and b must be shown together: with only FinnGen the DXA group looks
uniformly positive, and only HERMES reveals that it is not.

Run from the project root:  python analysis/scripts/make_paper_fig4.py
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

TAB = ROOT / "results" / "tables" / "15-modality"

US_KEYS = ["eBMD", "QUS_2", "QUS_3"]
DX_KEYS = ["FNBMD", "FNBMD_2", "FNBMD_3"]
SHORT = {"eBMD": "Morris 2018", "QUS_2": "PMID 34226706", "QUS_3": "PMID 30598549",
         "FNBMD": "PMID 39349818", "FNBMD_2": "same paper, 47k", "FNBMD_3": "PMID 41206123"}


SHORT_LABEL = {"eBMD": "Morris", "QUS_2": "34226706", "QUS_3": "30598549",
               "FNBMD": "39349818", "FNBMD_2": "47k", "FNBMD_3": "41206123"}


def forest(ax, d, pools, outcome, title, show_labels):
    d = d[d["outcome"] == outcome].set_index("exposure")
    us = [k for k in US_KEYS if k in d.index]
    dx = [k for k in DX_KEYS if k in d.index]
    rows = us + ["pool_us"] + dx + ["pool_dx"]

    y = np.arange(len(rows))[::-1]
    ax.axvline(0, color="0.45", lw=0.9, zorder=1)
    for yy, key in zip(y, rows):
        us_side = key in us or key == "pool_us"
        colour = pc.PROTECT if us_side else pc.WEAK
        if key.startswith("pool_"):
            grp = "ultrasound" if us_side else "dxa"
            p = pools[(pools["outcome"] == outcome) & (pools["group"] == grp)].iloc[0]
            ax.plot([p["lo"], p["hi"]], [yy, yy], color=colour, lw=1.9,
                    solid_capstyle="butt", zorder=2)
            ax.plot(p["b"], yy, "D", color="white", ms=5.8, mec=colour, mew=1.4, zorder=3)
            ax.annotate(f"{p['b']:+.3f}", (p["hi"], yy), textcoords="offset points",
                        xytext=(4, -1), va="center", fontsize=5.8, color=colour)
        else:
            s = d.loc[key]
            ax.plot([s["lo"], s["hi"]], [yy, yy], color=colour, lw=1.5,
                    solid_capstyle="butt", zorder=2)
            ax.plot(s["b"], yy, "o", color=colour, ms=4.8, mec="white", mew=0.55, zorder=3)
            ax.annotate(f"{s['b']:+.3f}", (s["hi"], yy), textcoords="offset points",
                        xytext=(4, -1), va="center", fontsize=5.8, color=colour)

    short_labels = []
    for k in rows:
        if k == "pool_us":
            short_labels.append(f"pooled, k={len(us)}")
        elif k == "pool_dx":
            short_labels.append(f"pooled, k={len(dx)}")
        else:
            tag = "US" if k in US_KEYS else "DXA"
            short_labels.append(f"{tag} {SHORT_LABEL.get(k, k)}")
    labels = []
    for k in rows:
        if k == "pool_us":
            labels.append(f"pooled, k={len(us)}")
        elif k == "pool_dx":
            labels.append(f"pooled, k={len(dx)}")
        elif k in us:
            labels.append(f"ultrasound, {SHORT[k]}")
        else:
            labels.append(f"DXA, {SHORT[k]}")
    ax.set_yticks(y)
    ax.set_yticklabels(labels if show_labels else short_labels, fontsize=6.0 if show_labels
                       else 5.6)
    for t, k in zip(ax.get_yticklabels(), rows):
        if k.startswith("pool_"):
            t.set_style("italic")
    ax.set_ylim(-1.6, len(rows) - 0.4)
    ax.margins(x=0.26)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.set_xlabel("effect on heart failure (log OR per SD)", fontsize=6.7)
    ax.set_title(title, loc="left", fontsize=8.0)

    c = pools[(pools["outcome"] == outcome)
              & (pools["group"] == "dxa_minus_ultrasound")].iloc[0]
    p_txt = "p<0.001" if c["p"] < 0.001 else f"p={c['p']:.2f}"
    ax.annotate(f"DXA - ultrasound\n{c['b']:+.3f} ({c['lo']:+.3f}, {c['hi']:+.3f}), {p_txt}",
                xy=(0.03, 0.03), xycoords="axes fraction", ha="left", va="bottom",
                fontsize=5.5, color="0.25")


def panel_c(ax, conc):
    d = conc.sort_values("b")
    y = np.arange(len(d))[::-1]
    ax.axvline(0, color="0.45", lw=0.9, zorder=1)
    ax.axvline(1, color="0.7", lw=0.9, ls=(0, (3, 2)), zorder=1)
    for yy, (_, s) in zip(y, d.iterrows()):
        within = (s["predictor"] in US_KEYS) == (s["predicted"] in US_KEYS)
        colour = pc.NULL if within else pc.POOL
        ax.plot([s["lo"], s["hi"]], [yy, yy], color=colour, lw=1.4,
                solid_capstyle="butt", zorder=2)
        ax.plot(s["b"], yy, "o", color=colour, ms=4.4, mec="white", mew=0.5, zorder=3)
        ax.annotate(f"{s['b']:+.2f}", (s["hi"], yy), textcoords="offset points",
                    xytext=(4, -1), va="center", fontsize=5.6, color=colour)
    ax.set_yticks(y)
    pretty = {"eBMD": "Morris", "QUS_2": "34226706", "QUS_3": "30598549",
              "FNBMD": "39349818", "FNBMD_2": "47k", "FNBMD_3": "41206123"}
    ax.set_yticklabels(
        [f'{pretty.get(r["predictor"], r["predictor"])} to '
         f'{pretty.get(r["predicted"], r["predicted"])}' for _, r in d.iterrows()],
        fontsize=5.6)
    ax.set_xlabel("SD of the predicted trait per SD of the predictor", fontsize=6.7)
    ax.set_ylim(-2.4, len(d) - 0.4)
    ax.margins(x=0.24)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_title("Both modalities measure the same trait", loc="left", fontsize=8.0)
    ax.annotate("grey: within a modality\nblack: ultrasound to DXA\ndashed: 1 SD per SD",
                xy=(0.03, 0.03), xycoords="axes fraction", ha="left", va="bottom",
                fontsize=5.4, color="0.3")


def main():
    df = pd.read_csv(TAB / "modality.csv")
    pools = pd.read_csv(TAB / "modality_pooled.csv")
    conc = pd.read_csv(TAB / "concordance.csv")

    pc.apply_figure_style(sizes=(8.0, 6.9, 5.9))
    fig = plt.figure(figsize=(7.4, 3.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 0.85, 1.0], wspace=1.05)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))

    forest(ax_a, df, pools, "FG_HF", "FinnGen, no UKB overlap", True)
    forest(ax_b, df, pools, "HF", "HERMES, shares UKB", False)
    panel_c(ax_c, conc)

    lo = min(ax.get_xlim()[0] for ax in (ax_a, ax_b))
    hi = max(ax.get_xlim()[1] for ax in (ax_a, ax_b))
    for ax in (ax_a, ax_b):
        ax.set_xlim(lo, hi)

    for ax, letter in zip([ax_a, ax_b, ax_c], "abc"):
        pc.panel_letter(ax, letter)

    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=pc.PROTECT, marker="o", ms=4.6, lw=1.5, label="heel ultrasound"),
               Line2D([], [], color=pc.WEAK, marker="o", ms=4.6, lw=1.5, label="DXA femoral neck"),
               Line2D([], [], color="white", marker="D", ms=5.0, mew=1.4, mec="0.35",
                      lw=1.6, label="pooled")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.075), ncol=3,
               frameon=False, fontsize=6.2, handlelength=1.7, handletextpad=0.55,
               columnspacing=1.8)

    pc.save(fig, "fig4_modality.png", ROOT)
    return 0


if __name__ == "__main__":
    sys.exit(main())

