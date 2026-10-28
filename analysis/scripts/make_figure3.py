"""Figure 3: why the transcriptomic layer cannot carry the comorbidity claim.

Panel a  DEG counts per dataset, split by direction, grouped by tissue side
Panel b  within-side replication: observed vs expected overlap for every dataset pair
Panel c  largest |log2FC| per dataset against its sample size
Panel d  Hallmark pathway counts and the bone-vascular convergence test

Run from the project root:  python analysis/scripts/make_figure3.py
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

DEG = ROOT / "results" / "tables" / "03-bulk-deg"
FIG = ROOT / "figures"

SIDE = {
    "GSE100927": "vascular", "GSE43292": "vascular", "GSE28829": "vascular",
    "GSE37558": "mechanism",
    "GSE56815": "bone", "GSE56814": "bone", "GSE7158": "bone", "GSE7429": "bone",
    "GSE13850": "bone", "GSE35956": "bone",
}
NSAMP = {"GSE100927": 104, "GSE43292": 64, "GSE28829": 29, "GSE56815": 80,
         "GSE56814": 73, "GSE7158": 26, "GSE7429": 20, "GSE13850": 40,
         "GSE35956": 10, "GSE37558": 32}
HUE = {"vascular": "#1f6fb4", "bone": "#b45f1f", "mechanism": "#3d8c6a"}


def load_deg():
    out = {}
    for f in sorted(DEG.glob("deg_GSE*.csv")):
        g = f.stem.replace("deg_", "")
        d = pd.read_csv(f)
        for c in ["logFC", "adj.P.Val"]:
            d[c] = pd.to_numeric(d[c], errors="coerce")
        out[g] = d
    return out


def panel_a(ax, deg):
    order = [g for g in deg if SIDE.get(g) == "vascular"] + \
            [g for g in deg if SIDE.get(g) == "mechanism"] + \
            [g for g in deg if SIDE.get(g) == "bone"]
    up = np.array([int(((deg[g]["adj.P.Val"] < 0.05) & (deg[g].logFC > 1)).sum()) for g in order])
    dn = np.array([int(((deg[g]["adj.P.Val"] < 0.05) & (deg[g].logFC < -1)).sum()) for g in order])
    y = np.arange(len(order))[::-1]
    ax.barh(y, up, height=0.62, color=[HUE[SIDE[g]] for g in order], edgecolor="none",
            label="up-regulated")
    ax.barh(y, -dn, height=0.62, color=[HUE[SIDE[g]] for g in order], edgecolor="none",
            alpha=0.42, label="down-regulated")
    ax.axvline(0, color="0.35", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=6.2)
    for lab, g in zip(ax.get_yticklabels(), order):
        lab.set_color(HUE[SIDE[g]])
    ax.set_xlabel("genes: FDR<0.05, |log2FC|>1", fontsize=6.9)
    lim = max(up.max(), dn.max()) * 1.18
    ax.set_xlim(-lim * 1.30, lim * 1.05)
    ax.margins(y=0.05)
    ax.legend(frameon=False, fontsize=6.0, loc="lower right", handlelength=1.1,
              borderpad=0.1, labelspacing=0.25)


def panel_b(ax):
    r = pd.read_csv(DEG / "within_side_replication.csv")
    r = r[(r.threshold == "fdr=0.05, lfc=0.5") & r.observed.notna() & r.expected.notna()]
    # a pair in which one dataset produced no DEGs at all says nothing about
    # reproducibility, so it is excluded rather than drawn as a failure
    r = r[r.expected >= 0.5].sort_values("fold")
    y = np.arange(len(r))
    ax.barh(y, r.fold, height=0.55, color=[HUE.get(s, "0.5") for s in r.side], edgecolor="none")
    ax.axvline(1, color="0.35", lw=0.9, ls=(0, (3, 2)))
    ax.set_yticks(y)
    ax.set_yticklabels([p.replace("_vs_", " vs ") for p in r.pair], fontsize=6.2)
    ax.set_xlabel("observed / expected DEG overlap", fontsize=6.9)
    ax.set_xscale("log")
    ax.set_xlim(0.5, 20)
    ax.set_ylim(-0.9, len(r) - 0.4)
    ax.annotate("chance", xy=(1, -0.75), ha="center", va="bottom", fontsize=5.9, color="0.35")
    ax.annotate("vascular comparisons\nclear chance; blood ones do not",
                xy=(0.97, 0.10), xycoords="axes fraction", ha="right", va="bottom",
                fontsize=6.0, color="0.25")


def panel_c(ax, deg):
    xs, ys, cs, labs = [], [], [], []
    for g, d in deg.items():
        if g in NSAMP:
            xs.append(NSAMP[g])
            ys.append(d.logFC.abs().max())
            cs.append(HUE[SIDE[g]])
            labs.append(g)
    ax.scatter(xs, ys, c=cs, s=30, edgecolor="white", linewidth=0.6, zorder=3)
    ax.text(0.015, 0.015, "blue vascular   orange bone/blood   green mechanism",
            transform=ax.transAxes, fontsize=5.8, color="0.30", va="bottom")
    offsets = {"GSE35956": (9, 3), "GSE56815": (6, 5), "GSE100927": (5, 6),
               "GSE43292": (7, -9), "GSE7158": (5, -8), "GSE56814": (-7, 6), "GSE28829": (7, -9),
               "GSE13850": (-8, -10), "GSE7429": (-7, -7), "GSE37558": (6, 2)}
    for x, y, l in zip(xs, ys, labs):
        ax.annotate(l, (x, y), textcoords="offset points", xytext=offsets.get(l, (5, 2)),
                    fontsize=5.6, color="0.30",
                    ha="right" if offsets.get(l, (5, 0))[0] < 0 else "left")
    ax.axhline(1.0, color="0.35", lw=0.9, ls=(0, (3, 2)))
    ax.annotate("2-fold", xy=(0.02, 1.04), xycoords=("axes fraction", "data"),
                ha="left", va="bottom", fontsize=5.9, color="0.35")
    ax.set_xscale("log")
    ax.set_xlabel("samples in the dataset (log)", fontsize=6.9)
    ax.set_ylabel("largest |log2FC| observed", fontsize=6.9)
    ax.margins(x=0.34, y=0.34)


def panel_d(ax):
    s = pd.read_csv(DEG / "pathway_summary.csv")
    agg = s.groupby(s["side"].astype(str))["n_pathways_fdr05"].sum()
    ax.bar(range(len(agg)), agg.values, width=0.55,
           color=[HUE.get(k, "0.5") for k in agg.index], edgecolor="none")
    ax.set_xticks(range(len(agg)))
    ax.set_xticklabels(agg.index, fontsize=6.8)
    ax.set_ylabel("enriched Hallmark\npathways (FDR<0.05)", fontsize=6.9)
    ax.set_ylim(0, max(agg.values) * 1.62)
    for i, v in enumerate(agg.values):
        ax.text(i, v + max(agg.values) * 0.04, str(int(v)), ha="center", fontsize=6.6,
                color="0.25")
    conv = pd.read_csv(DEG / "pathway_convergence.csv")
    bv = conv[conv.comparison == "bone vs vascular"].iloc[0]
    ax.annotate(f"shared = {int(bv.observed)}, expected {bv.expected:.1f}\npermutation p = {bv.perm_p:.2f}",
                xy=(0.5, 0.985), xycoords="axes fraction", ha="center", va="top",
                fontsize=6.1, color="#b03a2e")
    ax.margins(x=0.12)


def collide(fig):
    """Text overlaps, counting only artists actually drawn inside the canvas.

    Stale text objects left over from earlier axes still report a bbox (often 1x1 at the
    origin) and would otherwise be reported as overlaps.
    """
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
    apply_figure_style(sizes=(8.0, 6.9, 6.0))
    deg = load_deg()
    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.12, 1.0], hspace=0.52, wspace=0.52)
    ax_a, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    ax_c, ax_d = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    panel_a(ax_a, deg)
    panel_b(ax_b)
    panel_c(ax_c, deg)
    panel_d(ax_d)

    ax_a.set_title("Atherosclerosis moves hundreds of genes", loc="left", fontsize=8.0)
    ax_b.set_title("Only the vascular side reproduces", loc="left", fontsize=8.0)
    ax_c.set_title("Blood studies stay below 2-fold", loc="left", fontsize=8.0)
    ax_d.set_title("Pathways do not converge", loc="left", fontsize=8.0)
    for ax, letter in zip([ax_a, ax_b, ax_c, ax_d], "abcd"):
        panel_letter(ax, letter)

    for t in list(fig.findobj(matplotlib.text.Text)):
        if t.axes is None and t.get_figure() is fig:
            t.set_visible(False)
    FIG.mkdir(exist_ok=True)
    out = FIG / "fig3_deg_reproducibility.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"saved {out.relative_to(ROOT).as_posix()}")
    hits = collide(fig)
    print(f"collision scan: {len(hits)} overlaps")
    for h in hits[:12]:
        print("   ", ascii(h))
    return hits


if __name__ == "__main__":
    main()

