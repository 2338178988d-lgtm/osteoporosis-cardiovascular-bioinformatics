"""Figure 10: the exposure question 鈥?phenotype, or file?

Panel a  six bone-density files against FinnGen heart failure (no UKB overlap), grouped by
         modality, with the within-modality pool drawn as a diamond.
Panel b  the same six files against HERMES heart failure. Shown next to panel a on purpose:
         the DXA side is stable in neither direction across the two outcome cohorts, and a
         figure that showed only panel a would hide that.
Panel c  reciprocal prediction between the files. A file that does not predict the same
         phenotype measured elsewhere cannot be trusted in a and b, so this is the integrity
         check underneath the other two panels.

Run from the project root:  python analysis/scripts/make_figure10.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
from _paper_common import load_style
load_style()

TAB = ROOT / "results" / "tables" / "15-modality"
FIG = ROOT / "figures"

US = "#1f6fb4"
DXA = "#b45f1f"
GREY = "#8a8a8a"

SHORT = {
    "eBMD": ("heel US", "Morris 2018"),
    "QUS_2": ("heel US", "PMID 34226706"),
    "QUS_3": ("heel US", "PMID 30598549"),
    "FNBMD": ("DXA neck", "PMID 39349818"),
    "FNBMD_2": ("DXA neck", "same paper, 47k"),
    "FNBMD_3": ("DXA neck", "PMID 41206123"),
}
US_KEYS = ["eBMD", "QUS_2", "QUS_3"]
DX_KEYS = ["FNBMD", "FNBMD_2", "FNBMD_3"]


def forest(ax, df, pools, outcome, title, show_labels=True):
    d = df[df["outcome"] == outcome].set_index("exposure")
    us = [k for k in US_KEYS if k in d.index]
    dx = [k for k in DX_KEYS if k in d.index]
    rows = us + ["pool_ultrasound"] + dx + ["pool_dxa"]
    y = np.arange(len(rows))[::-1]

    ax.axvline(0, color="0.35", lw=0.9, zorder=1)
    for yy, key in zip(y, rows):
        us_side = key in us or key == "pool_ultrasound"
        colour = US if us_side else DXA
        if key.startswith("pool"):
            grp = "ultrasound" if us_side else "dxa"
            p = pools[(pools["outcome"] == outcome) & (pools["group"] == grp)].iloc[0]
            ax.plot([p["lo"], p["hi"]], [yy, yy], color=colour, lw=1.9,
                    solid_capstyle="butt", zorder=2)
            ax.plot(p["b"], yy, "D", color=GREY, ms=5.8, mec=colour, mew=1.3, zorder=3)
            ax.annotate(f"{p['b']:+.3f}", (p["hi"], yy), textcoords="offset points",
                        xytext=(4, -1), va="center", fontsize=5.8, color=colour)
        else:
            s = d.loc[key]
            ax.plot([s["lo"], s["hi"]], [yy, yy], color=colour, lw=1.5,
                    solid_capstyle="butt", zorder=2)
            ax.plot(s["b"], yy, "o", color=colour, ms=4.8, mec="white", mew=0.55, zorder=3)
            ax.annotate(f"{s['b']:+.3f}", (s["hi"], yy), textcoords="offset points",
                        xytext=(4, -1), va="center", fontsize=5.8, color=colour)

    labels = []
    for k in rows:
        if k.startswith("pool"):
            n = len(us) if k == "pool_ultrasound" else len(dx)
            labels.append(f"pooled k={n}")
        else:
            labels.append(f"{SHORT[k][0]}, {SHORT[k][1]}")
    ax.set_yticks(y)
    if show_labels:
        ax.set_yticklabels(labels, fontsize=6.0)
        for t, k in zip(ax.get_yticklabels(), rows):
            if k.startswith("pool"):
                t.set_style("italic")
                t.set_color(GREY)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel("effect on heart failure (log OR per SD)", fontsize=6.7)
    ax.set_ylim(-1.9, len(rows) - 0.4)
    ax.margins(x=0.22)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.set_title(title, loc="left", fontsize=7.6)

    c = pools[(pools["outcome"] == outcome)
              & (pools["group"] == "dxa_minus_ultrasound")].iloc[0]
    ax.annotate(f"DXA - US {c['b']:+.3f} ({c['lo']:+.3f}, {c['hi']:+.3f})\np = {c['p']:.1e}",
                xy=(0.03, 0.03), xycoords="axes fraction", ha="left", va="bottom",
                fontsize=5.6, color="0.25")


def concordance(ax, conc):
    d = conc.sort_values("b")
    y = np.arange(len(d))[::-1]
    ax.axvline(0, color="0.35", lw=0.9, zorder=1)
    ax.axvline(1, color="0.7", lw=0.8, ls=(0, (3, 2)), zorder=1)
    for yy, (_, s) in zip(y, d.iterrows()):
        colour = US if s["predicted"] in US_KEYS else DXA
        ax.plot([s["lo"], s["hi"]], [yy, yy], color=colour, lw=1.4,
                solid_capstyle="butt", zorder=2)
        ax.plot(s["b"], yy, "o", color=colour, ms=4.4, mec="white", mew=0.5, zorder=3)
        ax.annotate(f"{s['b']:+.2f}", (s["hi"], yy), textcoords="offset points",
                    xytext=(4, -1), va="center", fontsize=5.6, color=colour)
    ax.set_yticks(y)
    ax.set_yticklabels([r["pair_label"].replace("->", " to ") for _, r in d.iterrows()],
                       fontsize=5.8)
    ax.set_xlabel("SD of predicted trait per SD of predictor", fontsize=6.7)
    ax.set_ylim(-1.9, len(d) - 0.4)
    ax.margins(x=0.20)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.set_title("Same phenotype, same scale", loc="left", fontsize=7.6)
    ax.annotate("dashed line: one SD per SD (identical scale)",
                xy=(0.03, 0.10), xycoords="axes fraction", ha="left", va="bottom",
                fontsize=5.6, color="0.4")
    ax.annotate("coloured by the trait being predicted",
                xy=(0.03, 0.03), xycoords="axes fraction", ha="left", va="bottom",
                fontsize=5.6, color="0.25")


def main():
    df = pd.read_csv(TAB / "modality.csv")
    pools = pd.read_csv(TAB / "modality_pooled.csv")
    conc = pd.read_csv(TAB / "concordance.csv")

    apply_figure_style(sizes=(7.6, 6.6, 5.8))
    fig = plt.figure(figsize=(7.4, 3.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.80, 1.0], wspace=0.62)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))

    forest(ax_a, df, pools, "FG_HF", "FinnGen, no UKB overlap", show_labels=True)
    forest(ax_b, df, pools, "HF", "HERMES, shares UKB", show_labels=False)
    lo = min(ax_a.get_xlim()[0], ax_b.get_xlim()[0])
    hi = max(ax_a.get_xlim()[1], ax_b.get_xlim()[1])
    for ax in (ax_a, ax_b):
        ax.set_xlim(lo, hi)
    concordance(ax_c, conc)

    for ax, letter in zip([ax_a, ax_b, ax_c], "abc"):
        panel_letter(ax, letter)

    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=US, marker="o", ms=4.4, lw=1.5, label="heel ultrasound"),
               Line2D([], [], color=DXA, marker="o", ms=4.4, lw=1.5, label="DXA femoral neck"),
               Line2D([], [], color=GREY, marker="D", ms=5.0, lw=1.6, label="pooled")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.42, -0.10), ncol=3,
               frameon=False, fontsize=6.0, handlelength=1.7, handletextpad=0.55,
               columnspacing=1.5)

    FIG.mkdir(exist_ok=True)
    out = FIG / "fig10_modality.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"saved {out.relative_to(ROOT).as_posix()}")
    hits = collide(fig)
    print(f"collision scan: {len(hits)} overlaps")
    for h in hits[:10]:
        print("   ", ascii(h))
    return hits


def collide(fig):
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    box = fig.bbox
    texts = []
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_text().strip() or not t.get_visible():
            continue
        b = t.get_window_extent(rend)
        if b.width < 3 or b.height < 3:
            continue
        if b.x0 < box.x0 - 1 or b.x1 > box.x1 + 1 or b.y0 < box.y0 - 1 or b.y1 > box.y1 + 1:
            continue
        texts.append((t, b))
    ticks = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
             for ax in fig.axes}
    hits = []
    for i, (a, ba) in enumerate(texts):
        for b, bb in texts[i + 1:]:
            if ba.overlaps(bb):
                hits.append((a.get_text()[:30], b.get_text()[:30]))
    for t, bt in texts:
        if any(t in v for v in ticks.values()):
            continue
        for ax in fig.axes:
            for sp in ax.spines.values():
                if sp.get_visible() and bt.overlaps(sp.get_window_extent(rend)):
                    hits.append((t.get_text()[:30], "spine"))
    return hits


if __name__ == "__main__":
    main()
