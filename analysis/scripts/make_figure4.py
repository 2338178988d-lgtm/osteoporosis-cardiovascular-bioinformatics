"""Figure 4: stress-testing the BMD -> heart-failure signal (module 05-hf-followup).

Panel a  replication across three heart-failure GWAS
Panel b  multivariable MR: does the eBMD effect attenuate when a mediator is added?
Panel c  sensitivity: Steiger filtering, shared-locus removal, outlier removal

Run from the project root:  python analysis/scripts/make_figure4.py
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

HF = ROOT / "results" / "tables" / "05-hf-followup"
MR2 = ROOT / "results" / "tables" / "02-mr"
EXT = ROOT / "results" / "tables" / "08-external"
FIG = ROOT / "figures"

MAIN = "#1f6fb4"
ALT = "#b45f1f"
GREY = "#7a7a7a"


def panel_a(ax):
    d = pd.read_csv(HF / "fg_replication.csv")
    bbj = pd.read_csv(EXT / "ea_replication.csv").set_index("outcome").loc["FG_EA_HF"]
    d = pd.concat([d, bbj.to_frame().T.assign(outcome="FG_EA_HF")], ignore_index=True)
    order = ["HF", "HF_SHAH", "FG_HF", "FG_EA_HF"]
    label = {"HF": "HERMES\n47,309 cases",
             "HF_SHAH": "Shah 2020\n47,309 cases",
             "FG_HF": "FinnGen R12\n37,653 cases",
             "FG_EA_HF": "BioBank Japan\nchronic HF"}
    d = d.set_index("outcome").loc[order]
    y = np.arange(len(d))[::-1]
    ax.axvline(0, color="0.35", lw=0.8, zorder=1)
    for yy, (k, r) in zip(y, d.iterrows()):
        c = {"FG_HF": MAIN, "FG_EA_HF": ALT}.get(k, GREY)
        ax.plot([r["lo"], r["hi"]], [yy, yy], color=c, lw=1.7,
                alpha=1.0 if r["p"] < 0.05 else 0.62, solid_capstyle="butt", zorder=2)
        ax.plot(r["b"], yy, "o", color=c, ms=5.6, mec="white", mew=0.6,
                alpha=1.0 if r["p"] < 0.05 else 0.62, zorder=3)
        ax.annotate(f"OR {r['or_']:.3f}", (r["b"], yy), textcoords="offset points",
                    xytext=(0, 9), ha="center", fontsize=6.1, color=c)
    ax.set_yticks(y)
    ax.set_yticklabels([label[k] for k in d.index], fontsize=6.5)
    ax.set_xlabel("effect on heart failure (log OR per SD of BMD)", fontsize=6.9)
    ax.set_ylim(-1.75, len(d) - 0.55)
    ax.margins(x=0.22)
    ax.annotate("grey = shares UK Biobank with the exposure\nblue = no overlap, replicated\n"
                "orange = no overlap, 46 % power", xy=(0.02, 0.03), xycoords="axes fraction",
                ha="left", va="bottom", fontsize=5.8, color="0.30", linespacing=1.4)


def panel_b(ax):
    """FinnGen subtype test: present in all HF and obese HF, absent with coronary disease."""
    d = pd.read_csv(HF / "fg_subtype.csv")
    d = d.set_index("outcome").loc[["FG_HF", "FG_HF_OBESE", "FG_HF_CHD"]]
    label = {"FG_HF": "all heart failure\n37,653",
             "FG_HF_OBESE": "HF + BMI>=25\n22,965",
             "FG_HF_CHD": "HF + coronary disease\n17,532"}
    y = np.arange(len(d))[::-1]
    ax.axvline(0, color="0.35", lw=0.8, zorder=1)
    for yy, (k, r) in zip(y, d.iterrows()):
        sig = r["p"] < 0.05
        c = MAIN if sig else "#b03a2e"
        ax.plot([r["lo"], r["hi"]], [yy, yy], color=c, lw=1.7,
                alpha=1.0 if sig else 0.78, solid_capstyle="butt", zorder=2)
        ax.plot(r["b"], yy, "o", color=c, ms=5.6, mec="white", mew=0.6, zorder=3)
        ax.annotate(f"{r['b']:+.3f}", (r["b"], yy), textcoords="offset points",
                    xytext=(0, 9), ha="center", fontsize=6.1, color=c)
    ax.set_yticks(y)
    ax.set_yticklabels([label[k] for k in d.index], fontsize=6.5)
    ax.set_xlabel("effect on heart failure (log OR)", fontsize=6.9)
    ax.set_ylim(-0.85, len(d) - 0.35)
    ax.margins(x=0.16)
    ax.annotate("absent when coronary disease\nis part of the definition",
                xy=(0.02, 0.03), xycoords="axes fraction", ha="left", va="bottom",
                fontsize=5.9, color="#b03a2e", linespacing=1.4)

def panel_c(ax):
    rep = pd.read_csv(HF / "hf_replication.csv").set_index("outcome")
    sl = pd.read_csv(HF / "hf_shared_loci.csv").set_index("subset")
    sens = pd.read_csv(MR2 / "mr_pleiotropy_sensitivity.csv")
    s = sens[(sens.exposure == "eBMD") & (sens.outcome == "HF")].iloc[0]
    items = [
        ("all 592 instruments", float(rep.loc["HF", "b"]), float(rep.loc["HF", "se"]), GREY),
        ("outliers removed", float(s["trim_ivw_b"]), float(s["trim_ivw_se"]), MAIN),
        ("HF-near loci removed", float(sl.loc["excluding HF-near loci", "b"]),
         float(sl.loc["excluding HF-near loci", "se"]), MAIN),
    ]
    y = np.arange(len(items))[::-1]
    ax.axvline(0, color="0.35", lw=0.8, zorder=1)
    for yy, (lab, b, se, c) in zip(y, items):
        ax.plot([b - 1.96 * se, b + 1.96 * se], [yy, yy], color=c, lw=1.6,
                solid_capstyle="butt", zorder=2)
        ax.plot(b, yy, "o", color=c, ms=5.2, mec="white", mew=0.5, zorder=3)
        ax.annotate(f"{b:+.3f}", (b, yy), textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=6.2, color=c)
    ax.set_yticks(y)
    ax.set_yticklabels([i[0] for i in items], fontsize=6.6)
    ax.set_xlabel("eBMD effect on HF (log OR)", fontsize=6.9)
    ax.set_ylim(-0.75, len(items) - 0.35)
    ax.margins(x=0.10)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    ax.annotate("estimate strengthens under both filters", xy=(0.02, 0.03),
                xycoords="axes fraction", ha="left", va="bottom", fontsize=6.1, color="0.25")


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
    apply_figure_style(sizes=(8.0, 6.9, 6.0))
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.30, 1.10], wspace=0.72)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))

    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)

    ax_a.set_title("Third population inconclusive", loc="left", fontsize=7.6)
    ax_b.set_title("Absent from coronary heart failure", loc="left", fontsize=8.0)
    ax_c.set_title("Strengthens under filters", loc="left", fontsize=8.0)
    for ax, letter in zip([ax_a, ax_b, ax_c], "abc"):
        panel_letter(ax, letter)

    FIG.mkdir(exist_ok=True)
    out = FIG / "fig4_hf_followup.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"saved {out.relative_to(ROOT).as_posix()}")
    hits = collide(fig)
    print(f"collision scan: {len(hits)} overlaps")
    for h in hits[:10]:
        print("   ", ascii(h))
    return hits


if __name__ == "__main__":
    main()

