"""Figure 7: closing the two remaining alternatives 鈥?body composition, and pathway.

Panel a  eBMD effect on heart failure, unadjusted and adjusted for lean mass / height / BMI
Panel b  lean mass on its own instruments: marginal effect, then adjusted for BMI and height
Panel c  gene-set scan: observed largest |z| against the permutation null for both cohorts

Run from the project root:  python analysis/scripts/make_figure7.py
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

BC = ROOT / "results" / "tables" / "11-body-composition"
PS = ROOT / "results" / "tables" / "12-pathway-stratified"
FIG = ROOT / "figures"

MAIN = "#1f6fb4"
ALT = "#b45f1f"
GREY = "#8a8a8a"
BAD = "#b03a2e"


def panel_a(ax):
    d = pd.read_csv(BC / "bodycomp_mvmr.csv")
    labels = {"M1_eBMD": "eBMD alone", "M2_eBMD_ALM": "+ lean mass",
              "M3_eBMD_HEIGHT": "+ height", "M4_eBMD_ALM_BMI": "+ lean mass & BMI"}
    order = list(labels)
    y = np.arange(len(order))[::-1]
    ax.axvline(0, color="0.35", lw=0.9, zorder=1)
    for j, (coh, col) in enumerate([("HF", MAIN), ("FG_HF", ALT)]):
        for yy, m in zip(y, order):
            s = d[(d.outcome == coh) & (d.model == m)]
            if s.empty:
                continue
            s = s.iloc[0]
            off = 0.16 if j == 0 else -0.16
            ax.plot([s.direct_b - 1.96 * s.direct_se, s.direct_b + 1.96 * s.direct_se],
                    [yy + off] * 2, color=col, lw=1.4, solid_capstyle="butt", zorder=2)
            ax.plot(s.direct_b, yy + off, "o", color=col, ms=4.6, mec="white", mew=0.5, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([labels[m] for m in order], fontsize=6.8)
    ax.set_xlabel("eBMD effect on heart failure (log OR)", fontsize=6.9)
    ax.margins(x=0.18, y=0.10)
    ax.annotate("HERMES", xy=(0.02, 0.97), xycoords="axes fraction", ha="left", va="top",
                fontsize=6.2, color=MAIN, fontweight="bold")
    ax.annotate("FinnGen", xy=(0.02, 0.89), xycoords="axes fraction", ha="left", va="top",
                fontsize=6.2, color=ALT, fontweight="bold")
    ax.annotate("attenuation 鈮?9 % in all models", xy=(0.02, 0.02), xycoords="axes fraction",
                ha="left", va="bottom", fontsize=5.9, color="0.3")


def panel_b(ax):
    d = pd.read_csv(BC / "bodycomp_alm_own.csv")
    labels = {"univariable": "lean mass alone", "ALM+BMI": "+ BMI",
              "ALM+HEIGHT": "+ height", "ALM+BMI+HEIGHT": "+ BMI & height"}
    order = list(labels)
    y = np.arange(len(order))[::-1]
    ax.axvline(0, color="0.35", lw=0.9, zorder=1)
    for j, (coh, col) in enumerate([("HF", MAIN), ("FG_HF", ALT)]):
        for yy, m in zip(y, order):
            s = d[(d.outcome == coh) & (d.model == m)]
            if s.empty:
                continue
            s = s.iloc[0]
            off = 0.16 if j == 0 else -0.16
            sig = s.p < 0.05
            a = 1.0 if sig else 0.45
            ax.plot([s.b - 1.96 * s.se, s.b + 1.96 * s.se], [yy + off] * 2,
                    color=col, lw=1.4, alpha=a, solid_capstyle="butt", zorder=2)
            ax.plot(s.b, yy + off, "o", color=col, ms=4.6, mec="white", mew=0.5,
                    alpha=a, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([labels[m] for m in order], fontsize=6.8)
    ax.set_xlabel("lean-mass effect on heart failure (log OR)", fontsize=6.9)
    ax.margins(x=0.20, y=0.10)
    ax.annotate("faded = not significant", xy=(0.02, 0.02), xycoords="axes fraction",
                ha="left", va="bottom", fontsize=5.9, color="0.3")


def panel_c(ax):
    null = pd.read_csv(PS / "geneset_null.csv")
    fw = pd.read_csv(PS / "geneset_fwer.csv").set_index("outcome")
    bins = np.linspace(0, 8, 33)
    for coh, col in [("HF", MAIN), ("FG_HF", ALT)]:
        v = null[null.outcome == coh]["max_abs_z"]
        ax.hist(v, bins=bins, color=col, alpha=0.42, edgecolor="none", density=True)
        obs = float(fw.loc[coh, "observed_max_abs_z"])
        ax.axvline(obs, color=col, lw=1.4, ls="--", ymax=0.58, zorder=3)
    ax.set_xlabel("largest |z| across gene sets", fontsize=6.9)
    ax.set_ylabel("permutation density", fontsize=6.9)
    ax.margins(x=0.03, y=0.16)
    ax.annotate(f"observed HERMES {fw.loc['HF', 'observed_max_abs_z']:.2f}\n"
                f"observed FinnGen {fw.loc['FG_HF', 'observed_max_abs_z']:.2f}\n"
                f"family-wise p = {fw.loc['HF', 'p_familywise']:.2f} / "
                f"{fw.loc['FG_HF', 'p_familywise']:.2f}",
                xy=(0.97, 0.96), xycoords="axes fraction", ha="right", va="top",
                fontsize=5.9, color="0.2", linespacing=1.5)


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


def main():
    apply_figure_style(sizes=(8.0, 6.9, 5.9))
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.15], wspace=0.62)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))

    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)

    ax_a.set_title("Not muscle, not body size", loc="left", fontsize=8.0)
    ax_b.set_title("Lean mass is body size", loc="left", fontsize=8.0)
    ax_c.set_title("No gene set stands out", loc="left", fontsize=8.0)
    for ax, letter in zip([ax_a, ax_b, ax_c], "abc"):
        panel_letter(ax, letter)

    FIG.mkdir(exist_ok=True)
    out = FIG / "fig7_bodycomp_and_pathway.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"saved {out.relative_to(ROOT).as_posix()}")
    hits = collide(fig)
    print(f"collision scan: {len(hits)} overlaps")
    for h in hits[:10]:
        print("   ", ascii(h))
    return hits


if __name__ == "__main__":
    main()

