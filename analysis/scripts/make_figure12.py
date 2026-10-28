"""Figure 12: does the atrial-fibrillation and heart-failure effect survive other cohorts?

One claim: the effects are not a Finnish founder-population artefact and not an artefact of
sharing participants with the exposure 鈥?in particular, atrial fibrillation reproduces at
-0.071 (p = 0.0023) in a cohort that overlaps neither UK Biobank nor FinnGen.

Panel a  atrial fibrillation across five outcome cohorts
Panel b  heart failure across four outcome cohorts

Cohorts are drawn in one column, ordered by ancestry, with Finnish cohorts marked; the pooled
estimate is a diamond. The Finnish estimates are visibly the largest in both panels, which is
why the Finnish-versus-other contrast is printed rather than left to the eye.

Run from the project root:  python analysis/scripts/make_figure12.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
from _paper_common import load_style
load_style()

TAB = ROOT / "results" / "tables" / "14-cardiac-ladder"
FIG = ROOT / "figures"

FI = "#1f6fb4"        # Finnish
OTHER = "#8a8a8a"     # everything else
POOL = "#111111"


def pool(b, se):
    b, se = np.asarray(b, float), np.asarray(se, float)
    w = 1 / se ** 2
    bf = float(np.sum(w * b) / np.sum(w))
    q = float(np.sum(w * (b - bf) ** 2))
    dfree = len(b) - 1
    i2 = max(0.0, (q - dfree) / q) if q > 0 and dfree > 0 else 0.0
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (q - dfree) / c) if c > 0 and dfree > 0 else 0.0
    ws = 1 / (se ** 2 + tau2)
    br = float(np.sum(ws * b) / np.sum(ws))
    ser = float(np.sqrt(1 / np.sum(ws)))
    return br, ser, i2


SHORT = {
    "GCST90475970": "MVP",
    "GCST90624411": "EUR",
    "GCST90473571": "non-FI",
    "GCST90455657": "EUR",
    "GCST90473575": "non-FI",
    "GCST90475949": "MVP",
    "GCST90475951": "MVP-2",
    "GCST90473562": "non-FI ICD",
    "FinnGen strict": "FinnGen",
    "FinnGen reimbursement": "FinnGen",
    "HERMES": "HERMES",
    "FinnGen non-ischaemic": "FinnGen",
}


def short_label(label):
    """Compact row label: source and case count, nothing else."""
    n = label.split("(")[-1].rstrip(")") if "(" in label else ""
    for full, short in SHORT.items():
        if label.startswith(full):
            return f"{short} {n}" if n else short
    return label[:22]


def panel(ax, d, endpoint, title):
    """Every panel carries its own row labels: cohorts differ between endpoints, so a reader
    cannot map a grey row in one panel to a grey row in another."""
    d = d[d["endpoint"] == endpoint].sort_values("b")
    y = np.arange(len(d))[::-1]
    ax.axvline(0, color="0.45", lw=0.9, zorder=1)

    for yy, (_, s) in zip(y, d.iterrows()):
        fi = s["population"].startswith("Finnish")
        colour = FI if fi else OTHER
        lo, hi = s["b"] - 1.96 * s["se"], s["b"] + 1.96 * s["se"]
        ax.plot([lo, hi], [yy, yy], color=colour, lw=1.8 if fi else 1.4,
                solid_capstyle="butt", zorder=2)
        ax.plot(s["b"], yy, "o", color=colour, ms=5.2, mec="white", mew=0.6, zorder=3)
        ax.annotate(f"{s['b']:+.3f}", (hi, yy), textcoords="offset points",
                    xytext=(5, -1), va="center", fontsize=5.9, color=colour)

    b, se, i2 = pool(d["b"], d["se"])
    ypool = -1.0
    ax.plot([b - 1.96 * se, b + 1.96 * se], [ypool, ypool], color=POOL, lw=1.9,
            solid_capstyle="butt", zorder=2)
    ax.plot(b, ypool, "D", color=POOL, ms=5.8, mec="white", mew=0.7, zorder=3)
    ax.annotate(f"{b:+.3f}", (b + 1.96 * se, ypool), textcoords="offset points",
                xytext=(5, -1), va="center", fontsize=5.9, color=POOL)

    ys = list(y) + [ypool]
    ax.set_yticks(ys)
    ax.set_yticklabels([short_label(s["label"]) for _, s in d.iterrows()]
                       + [f"pooled, k={len(d)}"], fontsize=5.9)
    order = list(d.itertuples())
    for tick, yy in zip(ax.get_yticklabels(), ys):
        if yy == ypool:
            tick.set_style("italic")
            tick.set_color(POOL)
            continue
        idx = len(d) - 1 - int(yy)
        if 0 <= idx < len(order) and order[idx].population.startswith("Finnish"):
            tick.set_color(FI)

    ax.set_xlabel("effect of eBMD (log OR per SD)", fontsize=6.7)
    ax.set_ylim(-2.1, len(d) - 0.4)
    ax.margins(x=0.26)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.set_title(title, loc="left", fontsize=7.6)

    fin = d[d["population"].str.startswith("Finnish")]
    rest = d[~d["population"].str.startswith("Finnish")]
    if len(fin) and len(rest):
        bf, sef, _ = pool(fin["b"], fin["se"])
        br, ser, _ = pool(rest["b"], rest["se"])
        diff = bf - br
        sed = float(np.sqrt(sef ** 2 + ser ** 2))
        pv = float(2 * stats.norm.sf(abs(diff / sed)))
        p_txt = ("p < 0.001" if pv < 0.001 else
                 f"p = {pv:.3f}" if pv < 0.01 else f"p = {pv:.2f}")
        ax.annotate(f"Finnish {bf:+.3f} (k={len(fin)})\n"
                    f"other {br:+.3f} (k={len(rest)})\n"
                    f"difference {diff:+.3f}, {p_txt}",
                    xy=(0.02, 0.02), xycoords="axes fraction", ha="left", va="bottom",
                    fontsize=5.4, color="0.25")
    # the pooled-row statistics sit in the gap between the last cohort and the diamond
    ax.annotate(f"I$^2$ = {i2:.2f}", xy=(0.98, 0.215), xycoords="axes fraction",
                ha="right", va="center", fontsize=5.4, color="0.4")


def main():
    d = pd.read_csv(TAB / "big_cohorts_all.csv")

    apply_figure_style(sizes=(7.6, 6.6, 5.8))
    fig = plt.figure(figsize=(7.4, 4.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.0], wspace=0.90)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))

    panel(ax_a, d, "AF", "Atrial fibrillation")
    panel(ax_b, d, "HF", "Heart failure")
    panel(ax_c, d, "CM", "Non-ischaemic cardiomyopathy")

    lo = min(ax.get_xlim()[0] for ax in (ax_a, ax_b, ax_c))
    hi = max(ax.get_xlim()[1] for ax in (ax_a, ax_b, ax_c))
    for ax in (ax_a, ax_b, ax_c):
        ax.set_xlim(lo, hi)

    for ax, letter in zip([ax_a, ax_b, ax_c], "abc"):
        panel_letter(ax, letter)

    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=FI, marker="o", ms=4.8, lw=1.7, label="Finnish (FinnGen)"),
               Line2D([], [], color=OTHER, marker="o", ms=4.8, lw=1.4, label="other ancestry"),
               Line2D([], [], color=POOL, marker="D", ms=5.2, lw=1.8, label="pooled")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.075), ncol=3,
               frameon=False, fontsize=6.0, handlelength=1.7, handletextpad=0.55,
               columnspacing=1.6)

    FIG.mkdir(exist_ok=True)
    out = FIG / "fig12_cohort_replication.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"saved {out.relative_to(ROOT).as_posix()}")
    hits = collide(fig)
    print(f"collision scan: {len(hits)} overlaps")
    for h in hits[:8]:
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
