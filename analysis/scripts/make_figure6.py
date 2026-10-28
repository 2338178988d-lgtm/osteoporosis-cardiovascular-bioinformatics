"""Figure 6: from the genetic instrument to the cell that expresses it.

Panel a  the 589 instrument loci assigned to a nearest gene, and how many are measurable
Panel b  each instrument gene's preferred compartment, with the expression-matched null
Panel c  the split between the osteogenic transcription-factor programme and the
         bone-pathway genes, by cell type

Run from the project root:  python analysis/scripts/make_figure6.py
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

INT = ROOT / "results" / "tables" / "09-integration"
SCR = ROOT / "results" / "tables" / "07-scrna"
FIG = ROOT / "figures"

MAIN = "#1f6fb4"
ALT = "#b45f1f"
GREY = "#8a8a8a"
BAD = "#b03a2e"


def panel_a(ax):
    pref = pd.read_csv(INT / "gene_preferred_celltype.csv")
    n_total = len(pd.read_csv(INT / "instrument_genes.csv"))
    n_gene = len(pref)
    n_det = int((pref["preferred_celltype"] != "not_detected").sum())
    bars = [("with a nearest gene", n_gene, GREY),
            ("measurable in plaque", n_det, MAIN)]
    y = np.arange(len(bars))[::-1]
    for yy, (lab, v, c) in zip(y, bars):
        ax.barh(yy, v, height=0.5, color=c, edgecolor="none")
        ax.annotate(str(v), (v, yy), textcoords="offset points", xytext=(5, 0),
                    va="center", fontsize=7.0, color=c)
    ax.set_yticks(y)
    ax.set_yticklabels([b[0] for b in bars], fontsize=6.8)
    ax.set_xlim(0, n_total * 1.16)
    ax.set_xlabel(f"eBMD instruments (of {n_total})", fontsize=7.0)
    ax.margins(y=0.35)


def panel_b(ax):
    d = pd.read_csv(INT / "celltype_enrichment.csv").sort_values("fold")
    y = np.arange(len(d))
    for yy, (_, r) in zip(y, d.iterrows()):
        c = MAIN if r["fold"] > 1 and r["p"] < 0.05 else (BAD if r["fold"] < 0.9 else GREY)
        ax.barh(yy, r["fold"], height=0.58, color=c, edgecolor="none", alpha=0.90)
        ax.plot([r["null_q95"]], [yy], "|", color="0.20", ms=9, mew=1.3, zorder=4)
    ax.axvline(1, color="0.35", lw=0.9, ls=(0, (3, 2)))
    ax.set_yticks(y)
    ax.set_yticklabels([c.replace(" / ", "/") for c in d["celltype"]], fontsize=6.4)
    ax.set_xlabel("observed / expected share of instrument genes", fontsize=7.0)
    ax.set_xlim(0, 1.55)
    ax.margins(y=0.07)
    ax.annotate("tick = 95th percentile of the\nmatched null", xy=(0.99, 0.98),
                xycoords="axes fraction", ha="right", va="top", fontsize=5.9,
                color="0.30", linespacing=1.4)


def panel_c(ax):
    d = pd.read_csv(SCR / "programme_scores.csv")
    d = d.set_index("celltype")
    order = d.sort_values("osteogenic", ascending=False).index
    x = np.arange(len(order))
    w = 0.19
    for i, (col, lab, c) in enumerate([("osteogenic", "osteogenic TF", MAIN),
                                       ("calcification", "calcification", ALT),
                                       ("bmd_gene_score", "bone-pathway genes", "#3d8c6a")]):
        v = d.loc[order, col].to_numpy()
        ax.bar(x + (i - 1) * w, v, width=w, color=c, edgecolor="none", label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels([o.replace(" / ", "/") for o in order], rotation=90, ha="center",
                       fontsize=6.2)
    ax.set_ylabel("mean programme score", fontsize=7.0)
    ax.margins(y=0.16)
    ax.legend(frameon=False, fontsize=5.9, loc="upper right", handlelength=1.0,
              borderpad=0.1, labelspacing=0.25)


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
    fig = plt.figure(figsize=(7.2, 5.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.80, 1.05, 1.30], wspace=0.78)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))

    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)

    ax_a.set_title("Instruments trace to genes", loc="left", fontsize=8.0)
    ax_b.set_title("Those genes peak in fibroblasts", loc="left", fontsize=8.0)
    ax_c.set_title("Two bone-like programmes, two cells", loc="left", fontsize=8.0)
    for ax, letter in zip([ax_a, ax_b, ax_c], "abc"):
        panel_letter(ax, letter)

    FIG.mkdir(exist_ok=True)
    out = FIG / "fig6_instrument_to_cell.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"saved {out.relative_to(ROOT).as_posix()}")
    hits = collide(fig)
    print(f"collision scan: {len(hits)} overlaps")
    for h in hits[:10]:
        print("   ", ascii(h))
    return hits


if __name__ == "__main__":
    main()

