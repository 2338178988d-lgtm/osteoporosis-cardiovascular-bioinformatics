"""Figure 8: three robustness tests, and the one that changed the wording.

Panel a  the same effect estimated with two bone-density modalities - they disagree in sign
Panel b  six pre-specified negative-control endpoints for both modalities
Panel c  eBMD against left ventricular mass, in SD units

Run from the project root:  python analysis/scripts/make_figure8.py
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

RB = ROOT / "results" / "tables" / "13-robustness"
FIG = ROOT / "figures"

MAIN = "#1f6fb4"
ALT = "#b45f1f"
BAD = "#b03a2e"
GREY = "#8a8a8a"


def panel_a(ax):
    """Heel ultrasound vs DXA femoral neck, same outcome, same pipeline."""
    dxa = pd.read_csv(RB / "dxa_exposure.csv").set_index("outcome")
    items = [
        ("heel eBMD -> HERMES", -0.0336, 0.0130, MAIN),
        ("heel eBMD -> FinnGen", -0.0640, 0.0186, MAIN),
        ("DXA neck BMD -> HERMES", float(dxa.loc["HF", "b"]), float(dxa.loc["HF", "se"]), ALT),
        ("DXA neck BMD -> FinnGen", float(dxa.loc["FG_HF", "b"]),
         float(dxa.loc["FG_HF", "se"]), ALT),
    ]
    y = np.arange(len(items))[::-1]
    ax.axvline(0, color="0.35", lw=0.9, zorder=1)
    for yy, (lab, b, se, c) in zip(y, items):
        ax.plot([b - 1.96 * se, b + 1.96 * se], [yy, yy], color=c, lw=1.7,
                solid_capstyle="butt", zorder=2)
        ax.plot(b, yy, "o", color=c, ms=5.6, mec="white", mew=0.6, zorder=3)
        ax.annotate(f"{b:+.3f}", (b, yy), textcoords="offset points",
                    xytext=(-5 if b < 0 else 5, 9),
                    ha="right" if b < 0 else "left", fontsize=6.1, color=c)
    ax.set_yticks(y)
    ax.set_yticklabels([i[0] for i in items], fontsize=6.5)
    ax.set_xlabel("effect on heart failure (log OR)", fontsize=6.9)
    ax.set_ylim(-1.35, len(items) - 0.4)
    ax.margins(x=0.20)
    ax.annotate("same trait family, opposite sign", xy=(0.02, 0.03), xycoords="axes fraction",
                ha="left", va="bottom", fontsize=5.9, color=BAD)


def panel_b(ax):
    d = pd.read_csv(RB / "negative_controls.csv")
    neg = d[d.role == "negative_control"]
    labels = list(dict.fromkeys(neg.label))
    y = np.arange(len(labels))[::-1]
    ax.axvline(0, color="0.35", lw=0.9, zorder=1)
    for j, (exp, col) in enumerate([("eBMD", MAIN), ("FNBMD", ALT)]):
        for yy, lab in zip(y, labels):
            s = neg[(neg.exposure == exp) & (neg.label == lab)]
            if s.empty:
                continue
            s = s.iloc[0]
            off = 0.16 if j == 0 else -0.16
            sig = s.p < 0.05
            ax.plot([s.lo, s.hi], [yy + off] * 2, color=col, lw=1.3,
                    alpha=1.0 if sig else 0.55, solid_capstyle="butt", zorder=2)
            ax.plot(s.b, yy + off, "o", color=col, ms=4.2, mec="white", mew=0.4,
                    alpha=1.0 if sig else 0.55, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.set_xlabel("effect on negative control (log OR)", fontsize=6.9)
    ax.margins(x=0.16, y=0.10)
    ax.annotate("faded = p 鈮?0.05", xy=(0.02, 0.97), xycoords="axes fraction",
                ha="left", va="top", fontsize=5.9, color="0.3")
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=MAIN, marker="o", ms=4.0, lw=1.3, label="heel eBMD"),
               Line2D([], [], color=ALT, marker="o", ms=4.0, lw=1.3, label="DXA neck BMD")]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=6.0,
              handlelength=1.7, handletextpad=0.5, borderpad=0.15, labelspacing=0.35)


def panel_c(ax):
    d = pd.read_csv(RB / "lv_mass.csv")
    labels = {"LVM": "LV mass", "LVMI": "LV mass / BSA"}
    order = [k for k in ["LVM", "LVMI"] if k in set(d.outcome)]
    y = np.arange(len(order))[::-1]
    ax.axvline(0, color="0.35", lw=0.9, zorder=1)
    for yy, k in zip(y, order):
        s = d[d.outcome == k].iloc[0]
        sig = s.p < 0.05
        c = MAIN if sig else GREY
        ax.plot([s.lo_sd, s.hi_sd], [yy, yy], color=c, lw=1.7, solid_capstyle="butt",
                alpha=1.0 if sig else 0.7, zorder=2)
        ax.plot(s.b_sd, yy, "o", color=c, ms=5.6, mec="white", mew=0.6,
                alpha=1.0 if sig else 0.7, zorder=3)
        ax.annotate(f"{s.b_sd:+.3f}", (s.b_sd, yy), textcoords="offset points",
                    xytext=(0, 9), ha="center", fontsize=6.1, color=c)
    ax.set_yticks(y)
    ax.set_yticklabels([labels[k] for k in order], fontsize=6.5)
    ax.set_xlabel("LV mass per SD eBMD (SD)", fontsize=6.9)
    ax.set_ylim(-0.95, len(order) - 0.35)
    ax.margins(x=0.18)
    ax.annotate("591 instruments; power 66 %", xy=(0.02, 0.03),
                xycoords="axes fraction", ha="left", va="bottom", fontsize=5.9, color="0.3")


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
    fig = plt.figure(figsize=(7.8, 5.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.20, 1.20, 0.92], wspace=1.22)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))

    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)

    ax_a.set_title("Two modalities, opposite signs", loc="left", fontsize=8.0)
    ax_b.set_title("No generic pleiotropy", loc="left", fontsize=8.0)
    ax_c.set_title("Myocardium: small, negative", loc="left", fontsize=8.0)
    for ax, letter in zip([ax_a, ax_b, ax_c], "abc"):
        panel_letter(ax, letter)

    FIG.mkdir(exist_ok=True)
    out = FIG / "fig8_robustness.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"saved {out.relative_to(ROOT).as_posix()}")
    hits = collide(fig)
    print(f"collision scan: {len(hits)} overlaps")
    for h in hits[:10]:
        print("   ", ascii(h))
    return hits


if __name__ == "__main__":
    main()
