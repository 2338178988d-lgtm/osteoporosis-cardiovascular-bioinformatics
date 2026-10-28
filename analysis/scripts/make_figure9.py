"""Figure 9: the cardiac-load ladder 鈥?where on the path to heart failure the effect sits.

One claim: the bone-density effect is present on the myocardial rungs and absent on the load
and coronary rungs, so "the myocardium" is measured rather than inferred. The heart-failure
anchor is drawn as a band so each rung can be read against it.

Run from the project root:  python analysis/scripts/make_figure9.py
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

TAB = ROOT / "results" / "tables" / "14-cardiac-ladder"
FIG = ROOT / "figures"

PROTECT = "#1f6fb4"     # protective rung (negative estimate)
NULL = "#8a8a8a"        # compatible with no effect
ANCHOR = "#111111"      # the heart-failure anchor

# bottom of the ladder first; each entry is (outcome key, display label, group)
ORDER = [
    ("FG_HYPERTENS", "Hypertension", "load"),
    ("FG_MI_STRICT", "Myocardial infarction", "coronary"),
    ("FG_AF", "Atrial fibrillation", "rhythm"),
    ("FG_VHD", "Valvular heart disease", "valve"),
    ("FG_CAVS", "Calcific aortic stenosis", "valve"),
    ("FG_HYPERTENS_HD", "Hypertensive heart disease", "myocardium"),
    ("FG_NONISCH_CM", "Non-ischaemic cardiomyopathy", "myocardium"),
    ("FG_CARDMYO", "Cardiomyopathy (all)", "myocardium"),
    ("FG_HF", "Heart failure", "anchor"),
]


def main():
    d = pd.read_csv(TAB / "ladder.csv").set_index("outcome")
    apply_figure_style(sizes=(8.0, 6.9, 5.9))
    fig = plt.figure(figsize=(6.6, 4.5))
    ax = fig.add_subplot(111)

    rows = [r for r in ORDER if r[0] in d.index]
    y = np.arange(len(rows))[::-1]
    anchor = d.loc["FG_HF"]

    # the anchor CI as a band: every rung is read against the heart-failure estimate
    ax.axvspan(anchor["lo"], anchor["hi"], color=ANCHOR, alpha=0.07, lw=0, zorder=0)
    ax.axvline(0, color="0.35", lw=0.9, zorder=1)
    ax.axhline(anchor["b"] * 0 + (len(rows) - 1 - [r[0] for r in rows].index("FG_HF")) + 0.5,
               color="0.75", lw=0.7, zorder=1)

    for yy, (key, label, group) in zip(y, rows):
        s = d.loc[key]
        if key == "FG_HF":
            colour, ms, lw, hollow = ANCHOR, 6.4, 2.0, False
        elif s["p"] < 0.05 and s["b"] < 0:
            colour, ms, lw, hollow = PROTECT, 5.4, 1.8, False
        elif s["power_ref"] >= 0.65:
            colour, ms, lw, hollow = NULL, 4.8, 1.5, False
        else:
            # a null that could not have detected the heart-failure effect is not the same
            # statement as one that could: it is drawn open and labelled as such
            colour, ms, lw, hollow = NULL, 4.8, 1.3, True
        ax.plot([s["lo"], s["hi"]], [yy, yy], color=colour, lw=lw,
                solid_capstyle="butt", zorder=2,
                alpha=0.55 if hollow else 1.0)
        ax.plot(s["b"], yy, "o", color="white" if hollow else colour, ms=ms,
                mec=colour, mew=1.1 if hollow else 0.6, zorder=3)
        ax.annotate(f"{s['b']:+.3f}", (s["hi"], yy), textcoords="offset points",
                    xytext=(5, -1), va="center", fontsize=5.9, color=colour)

    ax.set_yticks(y)
    ax.set_yticklabels([r[1] for r in rows], fontsize=6.5)
    ax.set_xlabel("effect of eBMD on each endpoint (log OR per SD)", fontsize=6.9)
    ax.set_xlim(-0.16, 0.10)
    ax.set_ylim(-0.9, len(rows) - 0.35)
    ax.margins(x=0.02)

    # case counts, right-aligned in the margin so power is visible without another panel
    xmax = ax.get_xlim()[1]
    for yy, (key, _, _) in zip(y, rows):
        n = int(d.loc[key, "cases"])
        ax.annotate(f"{n / 1000:.0f}k", (xmax, yy), xycoords=("data", "data"),
                    textcoords="offset points", xytext=(30, -1), va="center",
                    ha="right", fontsize=5.6, color="0.45", annotation_clip=False)
    ax.annotate("cases", (xmax, len(rows) - 0.75), xycoords="data", textcoords="offset points",
                xytext=(30, 0), ha="right", fontsize=5.6, color="0.45",
                annotation_clip=False)

    ax.set_title("Null on blood pressure and coronary disease, protective signals across "
                 "non-coronary cardiac endpoints", loc="left", fontsize=8.0)
    panel_letter(ax, "a")

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [
        Line2D([], [], color=PROTECT, marker="o", ms=5.0, lw=1.8,
               label="protective, Bonferroni-significant"),
        Line2D([], [], color=NULL, marker="o", ms=4.6, lw=1.5, label="null, informative"),
        Line2D([], [], color=NULL, marker="o", ms=4.6, lw=1.3, mfc="white", mew=1.1,
               label="null, power < 65 %"),
        Line2D([], [], color=ANCHOR, marker="o", ms=5.8, lw=2.0, label="heart failure (anchor)"),
        Patch(facecolor=ANCHOR, alpha=0.10, edgecolor="none", label="heart-failure 95 % CI"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.45, -0.19), ncol=3,
              frameon=False, fontsize=6.0, handlelength=1.8, handletextpad=0.6,
              borderpad=0.1, labelspacing=0.35, columnspacing=1.3)

    FIG.mkdir(exist_ok=True)
    out = FIG / "fig9_cardiac_ladder.png"
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

