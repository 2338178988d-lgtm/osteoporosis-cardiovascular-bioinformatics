"""Figure 2 - atrial fibrillation and heart failure show concordant estimates across outcome datasets; the Finnish
point estimates are more negative in shared endpoint comparisons.

Panel a  atrial fibrillation across five outcome datasets or definitions, pooled, Finnish datasets marked
Panel b  heart failure across four outcome cohorts, pooled
Panel c  the Finnish-minus-other difference for the three endpoints measured in both settings

Claim: "Both effects are directionally concordant outside Finland; only atrial fibrillation has a cohort independent of both UK Biobank and FinnGen, and the Finnish estimate is more negative in every
shared endpoint comparison." Panel c reports the significance of those differences; it is not a
decoration and must not be dropped for space.

Run from the project root:  python analysis/scripts/make_paper_fig2.py
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import _paper_common as pc

pc.load_style()
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

TAB = ROOT / "results" / "tables" / "14-cardiac-ladder"

DIFF = "#2b6f6a"   # a difference between two pooled estimates

SHORT = {
    "GCST90475970": "MVP", "GCST90624411": "EUR", "GCST90473571": "non-Finnish",
    "GCST90455657": "EUR", "GCST90473575": "non-Finnish",
    "GCST90475949": "MVP-2", "GCST90475951": "MVP", "GCST90473562": "non-Finnish ICD",
    "FinnGen strict": "FinnGen", "FinnGen reimbursement": "FinnGen",
    "HERMES": "HERMES", "FinnGen non-ischaemic": "FinnGen",
}


def short(label):
    n = label.split("(")[-1].rstrip(")") if "(" in label else ""
    for full, s in SHORT.items():
        if label.startswith(full):
            return f"{s} {n}" if n else s
    return label[:24]


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


def cohort_panel(ax, d, endpoint, title):
    d = d[d["endpoint"] == endpoint].sort_values("b")
    y = np.arange(len(d))[::-1]
    ax.axvline(0, color="0.45", lw=0.9, zorder=1)

    lo_all = float((d["b"] - 1.96 * d["se"]).min())
    hi_all = float((d["b"] + 1.96 * d["se"]).max())
    pad = (hi_all - lo_all) * 0.22
    ax.set_xlim(lo_all - pad, hi_all + pad)
    for yy, (_, s) in zip(y, d.iterrows()):
        fi = s["population"].startswith("Finnish")
        colour = pc.FI if fi else pc.OTHER
        lo, hi = s["b"] - 1.96 * s["se"], s["b"] + 1.96 * s["se"]
        ax.plot([lo, hi], [yy, yy], color=colour, lw=1.8 if fi else 1.4,
                solid_capstyle="butt", zorder=2)
        ax.plot(s["b"], yy, "o", color=colour, ms=5.2, mec="white", mew=0.6, zorder=3)
        # put the value inside the panel when the interval reaches the right edge, otherwise it
        # spills into the next panel's tick labels
        if hi > ax.get_xlim()[1] * 0.55:
            ax.annotate(f"{s['b']:+.3f}", (lo, yy), textcoords="offset points",
                        xytext=(-5, -1), va="center", ha="right", fontsize=5.9, color=colour)
        else:
            ax.annotate(f"{s['b']:+.3f}", (hi, yy), textcoords="offset points",
                        xytext=(5, -1), va="center", fontsize=5.9, color=colour)

    b, se, i2 = pool(d["b"], d["se"])
    ypool = -1.0
    ax.plot([b - 1.96 * se, b + 1.96 * se], [ypool, ypool], color=pc.POOL, lw=1.9,
            solid_capstyle="butt", zorder=2)
    ax.plot(b, ypool, "D", color=pc.POOL, ms=5.8, mec="white", mew=0.7, zorder=3)
    ax.annotate(f"{b:+.3f}", (b + 1.96 * se, ypool), textcoords="offset points",
                xytext=(5, -1), va="center", fontsize=5.9, color=pc.POOL)

    ys = list(y) + [ypool]
    ax.set_yticks(ys)
    ax.set_yticklabels([short(s["label"]) for _, s in d.iterrows()]
                       + [f"pooled, k={len(d)}"], fontsize=6.0)
    order = list(d.itertuples())
    for tick, yy in zip(ax.get_yticklabels(), ys):
        if yy == ypool:
            tick.set_style("italic")
            tick.set_color(pc.POOL)
        else:
            idx = len(d) - 1 - int(yy)
            if 0 <= idx < len(order) and order[idx].population.startswith("Finnish"):
                tick.set_color(pc.FI)

    ax.set_xlabel("effect of eBMD (log OR per SD)", fontsize=6.7)
    ax.set_ylim(-2.0, len(d) - 0.4)
    ax.margins(x=0.26)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.set_title(title, loc="left", fontsize=8.0)
    ax.annotate(f"pooled I$^2$ = {i2:.2f}", xy=(0.02, 0.10), xycoords="axes fraction",
                ha="left", va="bottom", fontsize=5.5, color="0.4")


def contrast_panel(ax, allc):
    """Finnish minus non-Finnish, one row per endpoint.

    Every row is the same kind of quantity, so every row gets the same colour; significance is
    carried by the p-value printed beside the interval, not by hue. Using the blue that means
    "Finnish" elsewhere would make the one significant row look like the only Finnish one.
    """
    rows = []
    for endpoint, label in [("AF", "Atrial fibrillation"),
                            ("HF", "Heart failure"),
                            ("CM", "Non-ischaemic cardiomyopathy")]:
        d = allc[allc["endpoint"] == endpoint]
        fin, rest = d[d["finnish"]], d[~d["finnish"]]
        if not len(fin) or not len(rest):
            continue
        bf, sef, _ = pool(fin["b"], fin["se"])
        br, ser, _ = pool(rest["b"], rest["se"])
        diff = bf - br
        sed = float(np.sqrt(sef ** 2 + ser ** 2))
        rows.append(dict(label=label, bf=bf, br=br, diff=diff, sed=sed,
                         p=float(2 * stats.norm.sf(abs(diff / sed))),
                         k_fin=len(fin), k_other=len(rest)))

    y = np.arange(len(rows))[::-1]
    ax.axvline(0, color="0.45", lw=0.9, zorder=1)
    lo_all = min(r["diff"] - 1.96 * r["sed"] for r in rows)
    hi_all = max(r["diff"] + 1.96 * r["sed"] for r in rows)
    pad = (hi_all - lo_all) * 0.12
    ax.set_xlim(lo_all - pad, hi_all + (hi_all - lo_all) * 0.46)
    for yy, r in zip(y, rows):
        ax.plot([r["diff"] - 1.96 * r["sed"], r["diff"] + 1.96 * r["sed"]], [yy, yy],
                color=DIFF, lw=1.7, solid_capstyle="butt", zorder=2)
        ax.plot(r["diff"], yy, "s", color=DIFF, ms=4.8, mec="white", mew=0.6, zorder=3)
        p_txt = ("p<0.001" if r["p"] < 0.001
                 else f"p={r['p']:.3f}" if r["p"] < 0.01
                 else f"p={r['p']:.2f}")
        ax.annotate(f"{r['diff']:+.3f}, {p_txt}",
                    (r["diff"] + 1.96 * r["sed"], yy), textcoords="offset points",
                    xytext=(5, -1), va="center", fontsize=5.8, color=DIFF)
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=6.0)
    ax.set_xlabel("Finnish minus non-Finnish effect (log OR per SD)", fontsize=6.7)
    ax.set_ylim(-2.1, len(rows) - 0.4)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    ax.set_title("Finnish estimates are more negative;\n1 of 3 differences is significant", loc="left", fontsize=8.0)
    n_sig = sum(1 for r in rows if r["p"] < 0.05)
    ax.annotate(f"one of {len(rows)} reaches p < 0.05;\nthe other two are consistent with no "
                f"difference",
                xy=(0.02, 0.03), xycoords="axes fraction", ha="left", va="bottom",
                fontsize=5.4, color="0.4")


def main():
    allc = pd.read_csv(TAB / "big_cohorts_all.csv")

    pc.apply_figure_style(sizes=(8.0, 6.9, 5.9))
    fig = plt.figure(figsize=(7.8, 3.7))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.95, 1.0], wspace=1.70)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))

    cohort_panel(ax_a, allc, "AF", "Atrial fibrillation, 5 datasets/definitions")
    cohort_panel(ax_b, allc, "HF", "Heart failure, 4 cohorts")
    contrast_panel(ax_c, allc)

    lo = min(ax.get_xlim()[0] for ax in (ax_a, ax_b))
    hi = max(ax.get_xlim()[1] for ax in (ax_a, ax_b))
    for ax in (ax_a, ax_b):
        ax.set_xlim(lo, hi)

    for ax, letter in zip([ax_a, ax_b, ax_c], "abc"):
        pc.panel_letter(ax, letter)

    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=pc.FI, marker="o", ms=4.8, lw=1.8, label="Finnish (FinnGen)"),
               Line2D([], [], color=pc.OTHER, marker="o", ms=4.8, lw=1.4, label="non-Finnish cohorts"),
               Line2D([], [], color=pc.POOL, marker="D", ms=5.2, lw=1.9, label="pooled"),
               Line2D([], [], color=DIFF, marker="s", ms=4.8, lw=1.7,
                      label="Finnish - non-Finnish (panel c)")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.095), ncol=4,
               frameon=False, fontsize=6.2, handlelength=1.7, handletextpad=0.55,
               columnspacing=1.8)

    pc.save(fig, "fig2_replication.png", ROOT)
    return 0


if __name__ == "__main__":
    sys.exit(main())

