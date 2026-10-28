"""Figure 2: genetic sharing between bone mineral density / osteoporosis and
cardiovascular endpoints (module 02-mr), plus the transcriptomic non-convergence
(module 03-bulk-deg).

Panel a  forest of the two-sample MR estimates (LD-clumped instruments, IVW random)
Panel b  the four estimators for the three eBMD pairs, testing method dependence
Panel c  instrument-set accounting: candidates -> LD-independent -> used
Panel d  two-tissue DEG intersection vs chance expectation

Run from the project root with the project .venv:
    python analysis/scripts/make_figure2.py
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

MR = ROOT / "results" / "tables" / "02-mr"
DEG = ROOT / "results" / "tables" / "03-bulk-deg"
FIG = ROOT / "figures"

SHORT = {"eBMD": "BMD", "OP": "OP", "CAD": "CAD", "MI": "MI", "STROKE": "Stroke", "HF": "HF"}
ORDER = [("eBMD", "CAD"), ("eBMD", "MI"), ("eBMD", "STROKE"), ("eBMD", "HF"),
         ("eBMD", "OP"), ("OP", "CAD"), ("OP", "MI"), ("OP", "STROKE"), ("OP", "HF"),
         ("CAD", "eBMD"), ("MI", "eBMD"), ("STROKE", "eBMD"), ("HF", "eBMD"),
         ("OP", "eBMD"), ("CAD", "OP"), ("MI", "OP"), ("STROKE", "OP"), ("HF", "OP")]
BONE = {"eBMD", "OP"}


def load_mr():
    d = pd.concat([pd.read_csv(MR / "mr_results_forward.csv"),
                   pd.read_csv(MR / "mr_results_reverse.csv")], ignore_index=True)
    for c in ["b", "se", "p", "q_p", "intercept_p", "nsnp"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d


def forest(d):
    rows = []
    for e, o in ORDER:
        s = d[(d.exposure == e) & (d.outcome == o) & (d.method == "IVW (random)")]
        if not s.empty:
            s = s.iloc[0]
            rows.append(dict(exposure=e, outcome=o, b=s["b"], se=s["se"], p=s["p"], nsnp=s["nsnp"]))
    return pd.DataFrame(rows)


def panel_a(ax, r, sens):
    """Forest plot. Filled = survives pleiotropy correction; a cross means the nominally
    significant IVW estimate lost significance once outlier variants were removed, so the
    panel cannot show an effect the text rejects."""
    y = np.arange(len(r))[::-1]
    to_cvd = r.exposure.isin(BONE) & ~r.outcome.isin(BONE)
    to_bone = r.outcome.isin(BONE) & ~r.exposure.isin(BONE)
    col = np.where(to_cvd, "#1f6fb4", np.where(to_bone, "#b45f1f", "#5c5c5c"))
    ax.axvline(0, color="0.35", lw=0.8, zorder=1)
    for i, (yy, row) in enumerate(zip(y, r.itertuples())):
        s = sens.get((row.exposure, row.outcome), {})
        trim_p = s.get("trim_ivw_p", np.nan)
        nominal = row.p < 0.05
        robust = nominal and (not np.isfinite(trim_p) or trim_p < 0.05)
        weak = not np.isfinite(trim_p) and nominal      # estimate undefined after trimming
        alpha = 1.0 if robust else 0.38
        ax.plot([row.b - 1.96 * row.se, row.b + 1.96 * row.se], [yy, yy],
                color=col[i], lw=1.5, alpha=alpha, solid_capstyle="butt", zorder=2)
        ax.plot(row.b, yy, "o" if row.nsnp >= 30 else "D", color=col[i],
                ms=5.0 if robust else 4.0, mec="white", mew=0.5, alpha=alpha, zorder=3)
        if nominal and not robust:
            ax.plot(row.b, yy, "x", color="#b03a2e", ms=5.6, mew=1.1, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{SHORT[t.exposure]} 鈫?{SHORT[t.outcome]}" for t in r.itertuples()],
                       fontsize=6.9)
    ax.set_xlabel("MR effect (log OR, 95% CI)", fontsize=7.0, labelpad=4)
    ax.margins(x=0.08, y=0.02)
    for lab, grp in zip(ax.get_yticklabels(), np.where(to_cvd, "b", np.where(to_bone, "o", "g"))):
        lab.set_color({"b": "#1f6fb4", "o": "#b45f1f", "g": "#5c5c5c"}[grp])
    ax.annotate("blue bone鈫抍ardio 路 orange cardio鈫抌one 路 grey control\n"
                "filled = survives pleiotropy correction;  red 脳 = loses it;  鈼?= <30 instruments",
                xy=(1.0, -0.07), xycoords="axes fraction", ha="right", va="top",
                fontsize=5.8, color="0.30", linespacing=1.35)


def panel_b(ax, d):
    pairs = [("eBMD", "HF"), ("eBMD", "CAD"), ("eBMD", "OP")]
    methods = ["IVW (random)", "MR-Egger", "Weighted median", "Weighted mode"]
    short = ["IVW", "Egger", "W. median", "W. mode"]
    hues = ["#1f6fb4", "#8c8c8c", "#b45f1f", "#3d8c6a"]
    h = 0.19
    for i, m in enumerate(methods):
        for j, (e, o) in enumerate(pairs):
            s = d[(d.exposure == e) & (d.outcome == o) & (d.method == m)]
            if s.empty:
                continue
            s = s.iloc[0]
            yy = (len(pairs) - 1 - j) + (1.5 - i) * h
            ax.plot([s["b"] - 1.96 * s["se"], s["b"] + 1.96 * s["se"]], [yy, yy],
                    color=hues[i], lw=1.2, solid_capstyle="butt")
            ax.plot(s["b"], yy, "o", color=hues[i], ms=4.2, mec="white", mew=0.4,
                    label=short[i] if j == 0 else None)
    ax.axvline(0, color="0.35", lw=0.8, zorder=0)
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels([f"BMD 鈫?{SHORT[o]}" for _, o in pairs][::-1], fontsize=7.0)
    ax.set_xlabel("MR effect (log OR, 95% CI)", fontsize=7.0)
    ax.set_ylim(-0.55, len(pairs) - 0.10)
    ax.margins(x=0.08)
    ax.legend(frameon=False, fontsize=6.4, ncol=2, loc="lower left",
              handletextpad=0.3, columnspacing=0.9, borderpad=0.1)


def panel_c(ax):
    s = pd.read_csv(MR / "ld_clump_summary.csv").set_index("dataset")
    s = s.loc[["eBMD", "CAD", "MI", "STROKE", "HF", "OP"]].reset_index()
    y = np.arange(len(s))[::-1]
    ax.barh(y, s.candidates, height=0.50, color="0.83", edgecolor="none",
            label="p<5e-8 candidates")
    ax.barh(y, s.clumped, height=0.50, color="#1f6fb4", edgecolor="none",
            label="LD-independent")
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels([SHORT[k] for k in s.dataset], fontsize=6.9)
    ax.set_xlabel("variants (log scale)", fontsize=7.0)
    for yy, row in zip(y, s.itertuples()):
        ax.text(row.clumped * 1.35, yy, str(row.clumped), va="center", fontsize=6.3, color="#1f6fb4")
    ax.set_xlim(10, s.candidates.max() * 4)
    ax.set_ylim(-0.55, len(s) + 0.45)
    ax.margins(y=0.05)
    ax.legend(frameon=False, fontsize=6.1, loc="lower left", handlelength=1.1, borderpad=0.1)


def panel_d(ax):
    g = pd.read_csv(DEG / "intersection_grid.csv")
    g = g[(g.model == "bone_disc_vs_vascular_disc") & g.observed.notna() & g.expected.notna() & (g.n_b > 0)]
    labels = [t.replace("fdr=", "FDR<").replace("lfc=", "FC>").replace("0.0", "0")
              .replace(", ", "\n") for t in g.threshold]
    x = np.arange(len(g))
    ax.bar(x - 0.19, g.expected, width=0.36, color="0.82", edgecolor="none", label="chance")
    ax.bar(x + 0.19, g.observed, width=0.36, color="#1f6fb4", edgecolor="none", label="observed")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=0, fontsize=5.6)
    ax.set_ylabel("shared genes", fontsize=7.0)
    ax.set_xlabel("selection threshold", fontsize=7.0, labelpad=6)
    ax.margins(y=0.16)
    ax.legend(frameon=False, fontsize=6.2, loc="upper right", handlelength=1.1, borderpad=0.1)
    ax.annotate("below chance at every threshold", xy=(0.02, 0.97), xycoords="axes fraction",
                ha="left", va="top", fontsize=6.3, color="#b03a2e")


def collide(fig):
    """Text-vs-text and text-vs-foreign-spine overlaps; own tick labels are exempt."""
    rend = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(rend)) for t in fig.findobj(matplotlib.text.Text)
             if t.get_text().strip() and t.get_visible()]
    ticks = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
             for ax in fig.axes}
    hits = []
    for i, (a, ba) in enumerate(texts):
        for b, bb in texts[i + 1:]:
            if ba.overlaps(bb):
                hits.append((a.get_text()[:34], b.get_text()[:34]))
    for t, bt in texts:
        if any(t in v for v in ticks.values()):
            continue
        for ax in fig.axes:
            for s in ax.spines.values():
                if s.get_visible() and bt.overlaps(s.get_window_extent(rend)):
                    hits.append((t.get_text()[:34], "spine"))
    return hits


def main():
    apply_figure_style(sizes=(8.0, 6.9, 6.0))
    d = load_mr()
    fig = plt.figure(figsize=(7.2, 7.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1.0], height_ratios=[1.55, 1.0],
                          hspace=0.40, wspace=0.40)
    ax_a, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    ax_c, ax_d = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    r = forest(d)
    sens_df = pd.read_csv(MR / "mr_pleiotropy_sensitivity.csv")
    sens = {(row.exposure, row.outcome): dict(trim_ivw_p=row.trim_ivw_p,
                                              n_outliers=row.n_outliers)
            for row in sens_df.itertuples()}
    panel_a(ax_a, r, sens)
    panel_b(ax_b, d)
    panel_c(ax_c)
    panel_d(ax_d)

    ax_a.set_title("Only heart failure survives pleiotropy correction", loc="left", fontsize=8.0)
    ax_b.set_title("Direction holds across estimators", loc="left", fontsize=8.0)
    ax_c.set_title("Instruments after LD clumping", loc="left", fontsize=8.0)
    ax_d.set_title("Bone and vascular DEGs: no convergence", loc="left", fontsize=8.0)
    for ax, letter in zip([ax_a, ax_b, ax_c, ax_d], "abcd"):
        panel_letter(ax, letter)

    FIG.mkdir(exist_ok=True)
    out = FIG / "fig2_mr_and_deg.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"saved {out.relative_to(ROOT).as_posix()}")
    hits = collide(fig)
    print(f"collision scan: {len(hits)} overlaps")
    for h in hits[:12]:
        print("   ", ascii(h))
    save_panel_crops(str(out), panel_crops(fig))
    return hits


if __name__ == "__main__":
    main()

