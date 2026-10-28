"""Figure 5: ten circulating candidates, none of which carries the effect.

Panel a  attenuation of the eBMD -> HF effect by each mediator, both cohorts
Panel b  each mediator's own effect on HF, on its own instruments, both cohorts
Panel c  colocalisation posterior at the 600 instrument loci

Estradiol is drawn in panel a only: 13 autosomal instruments give its own-effect
estimate an SE of 0.36, which would set the scale of panel b while carrying no
information.

Run from the project root:  python analysis/scripts/make_figure5.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
from _paper_common import load_style
load_style()

AX = ROOT / "results" / "tables" / "06-bone-axis"
COL = ROOT / "results" / "tables" / "09-integration"
FIG = ROOT / "figures"

MAIN = "#1f6fb4"
ALT = "#b45f1f"
GREY = "#7d7d7d"

ORDER = ["FGF23", "SOST", "BGLAP", "DKK1", "PTH", "TESTO", "ESTRADIOL", "VITD", "CST3", "CRP"]
LABEL = {"FGF23": "FGF23", "SOST": "sclerostin", "BGLAP": "osteocalcin", "DKK1": "DKK1",
         "PTH": "PTH", "TESTO": "testosterone", "ESTRADIOL": "estradiol",
         "VITD": "25(OH)D", "CST3": "cystatin C", "CRP": "CRP"}
PANEL_B_SKIP = {"ESTRADIOL"}


def panel_a(ax):
    d = pd.read_csv(AX / "axis_mvmr.csv")
    piv = d.pivot_table(index="mediator", columns="outcome", values="attenuation_pct")
    y = np.arange(len(ORDER))[::-1]
    ax.axvline(0, color="0.35", lw=0.9, zorder=1)
    for j, (coh, col) in enumerate([("HF", MAIN), ("FG_HF", ALT)]):
        for yy, m in zip(y, ORDER):
            if m not in piv.index or coh not in piv.columns:
                continue
            v = piv.loc[m, coh]
            ax.plot(v, yy + (0.16 if j == 0 else -0.16), "o", color=col, ms=5.0,
                    mec="white", mew=0.5, zorder=3)
    ax.axvspan(-10, 10, color="0.90", zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[m] for m in ORDER], fontsize=7.0)
    ax.set_xlabel("share of the eBMD effect absorbed (%)", fontsize=7.0)
    ax.set_xlim(-16, 16)
    ax.margins(y=0.06)
    ax.annotate("grey band = 卤10%", xy=(0.98, 0.03), xycoords="axes fraction",
                ha="right", va="bottom", fontsize=6.0, color="0.35")
    ax.annotate("HERMES", xy=(0.02, 0.98), xycoords="axes fraction", ha="left", va="top",
                fontsize=6.3, color=MAIN, fontweight="bold")
    ax.annotate("FinnGen", xy=(0.02, 0.91), xycoords="axes fraction", ha="left", va="top",
                fontsize=6.3, color=ALT, fontweight="bold")


def panel_b(ax):
    d = pd.read_csv(AX / "axis_univariable.csv")
    y = np.arange(len(ORDER))[::-1]
    ax.axvline(0, color="0.35", lw=0.9, zorder=1)
    for j, (coh, col) in enumerate([("HF", MAIN), ("FG_HF", ALT)]):
        for yy, m in zip(y, ORDER):
            if m in PANEL_B_SKIP:
                continue
            s = d[(d.mediator == m) & (d.outcome == coh)]
            if s.empty:
                continue
            s = s.iloc[0]
            ax.plot([s["lo"], s["hi"]], [yy + (0.16 if j == 0 else -0.16)] * 2,
                    color=col, lw=1.3, alpha=0.85, solid_capstyle="butt", zorder=2)
            ax.plot(s["b"], yy + (0.16 if j == 0 else -0.16), "o", color=col, ms=4.4,
                    mec="white", mew=0.4, alpha=0.95, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[m] for m in ORDER], fontsize=7.0)
    ax.set_xlabel("mediator effect on heart failure (log OR)", fontsize=7.0)
    ax.annotate("estradiol omitted\n(SE 0.36, 13 SNPs)", xy=(0.02, 0.03),
                xycoords="axes fraction", ha="left", va="bottom", fontsize=5.8, color="0.35")
    ax.margins(x=0.16, y=0.06)


def panel_c(ax):
    d = pd.read_csv(COL / "coloc_loci.csv").dropna(subset=["PP_H4"])
    ax.hist(d["PP_H4"], bins=np.linspace(0, 1, 41), color=MAIN, edgecolor="none")
    ax.set_yscale("log")
    ax.set_xlabel("posterior probability of a shared causal variant (PP.H4)", fontsize=7.0)
    ax.set_ylabel("instrument loci (log scale)", fontsize=7.0)
    n4 = int((d["PP_H4"] >= 0.8).sum())
    med = float(d["PP_H4"].median())
    ax.annotate(f"median PP.H4 = {med:.3f}\n{n4} of {len(d)} loci 鈮?0.8",
                xy=(0.97, 0.96), xycoords="axes fraction", ha="right", va="top",
                fontsize=6.5, color="0.20", linespacing=1.5)
    ax.annotate("the one locus is ABO", xy=(0.98, 0.62), xycoords="axes fraction",
                ha="right", va="top", fontsize=6.2, color="#b03a2e")
    ax.margins(x=0.03, y=0.14)


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
    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.35, 1.0], wspace=0.62)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))

    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)

    ax_a.set_title("No mediator absorbs the effect", loc="left", fontsize=8.0)
    ax_b.set_title("None moves heart failure consistently", loc="left", fontsize=8.0)
    ax_c.set_title("Almost no locus colocalises", loc="left", fontsize=8.0)
    for ax, letter in zip([ax_a, ax_b, ax_c], "abc"):
        panel_letter(ax, letter)

    FIG.mkdir(exist_ok=True)
    out = FIG / "fig5_mediators_and_coloc.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"saved {out.relative_to(ROOT).as_posix()}")
    hits = collide(fig)
    print(f"collision scan: {len(hits)} overlaps")
    for h in hits[:10]:
        print("   ", ascii(h))
    return hits


if __name__ == "__main__":
    main()

