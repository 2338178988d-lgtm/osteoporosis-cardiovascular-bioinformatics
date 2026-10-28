"""Figure 11: what the extended panel decides about the ladder.

Four claims, one per panel, all from the same 600 instruments and the same Finnish cohort:

a  The effect is cardiac, not systemic matrix. Four venous/aortic endpoints are null while the
   three cardiac anchors are not.
b  Within the heart it is atrial fibrillation, not the conduction system and not arrhythmia in
   general 鈥?including the endpoint with 95 % power to have seen it.
c  The effect grows as the endpoint becomes more specific, which is what the dilution
   synthesis requires: no effect in all-cardiovascular or ischaemic disease, the full effect in
   non-ischaemic cardiomyopathy.
d  Valve disease: the ICD aggregate is protective, the surgical calcific subset and a direct
   calcification GWAS are not significant at 47-53 % power, so the calcific reading stays open.

Run from the project root:  python analysis/scripts/make_figure11.py
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

PROTECT = "#1f6fb4"
NULL = "#8a8a8a"
OPEN = "#b45f1f"


def rows_for(keys, labels, panel):
    d = panel.set_index("outcome")
    lad = pd.read_csv(TAB / "ladder.csv").set_index("outcome")
    out = []
    for k, lab in zip(keys, labels):
        src = d if k in d.index else lad
        if k not in src.index:
            continue
        r = src.loc[k]
        out.append(dict(key=k, label=lab, b=float(r["b"]), lo=float(r["lo"]),
                        hi=float(r["hi"]), p=float(r["p"]),
                        power=float(r["power_ref"]), cases=int(r["cases"])))
    return out


def forest(ax, items, rules=None):
    """One row per estimate, significance driving the mark, power printed in the margin."""
    y = np.arange(len(items))[::-1]
    ax.axvline(0, color="0.35", lw=0.9, zorder=1)
    for yy, s in zip(y, items):
        sig = s["p"] < 0.05
        if not sig and s["power"] < 0.65:
            colour, filled = OPEN, False
        elif sig and s["b"] < 0:
            colour, filled = PROTECT, True
        else:
            colour, filled = NULL, True
        ax.plot([s["lo"], s["hi"]], [yy, yy], color=colour, lw=1.5 if filled else 1.2,
                solid_capstyle="butt", zorder=2, alpha=1.0 if filled else 0.6)
        ax.plot(s["b"], yy, "o", color="white" if not filled else colour, ms=4.8,
                mec=colour, mew=1.2 if not filled else 0.55, zorder=3)
        ax.annotate(f"{s['b']:+.3f}", (s["hi"], yy), textcoords="offset points",
                    xytext=(4, -1), va="center", fontsize=5.6, color=colour)
        ax.annotate(f"{s['power']:.0%}", (1.0, yy), xycoords=("axes fraction", "data"),
                    textcoords="offset points", xytext=(16, -1), va="center",
                    ha="right", fontsize=5.3, color="0.5", annotation_clip=False)
    ax.set_yticks(y)
    ax.set_yticklabels([s["label"] for s in items], fontsize=6.0)
    ax.set_ylim(-1.5, len(items) - 0.4)
    ax.margins(x=0.20)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    if rules:
        for ypos, text in rules:
            ax.axhline(ypos, color="0.85", lw=0.7, zorder=0)
            ax.annotate(text, (0.01, ypos), xycoords=("axes fraction", "data"),
                        textcoords="offset points", xytext=(0, 3), ha="left", va="bottom",
                        fontsize=5.4, color="0.45", style="italic")
    ax.annotate("power", (1.0, len(items) - 0.7), xycoords=("axes fraction", "data"),
                textcoords="offset points", xytext=(16, 0), ha="right", fontsize=5.3,
                color="0.5", annotation_clip=False)


def panel_a(ax, panel):
    items = rows_for(
        ["FG_VARICVE", "FG_VTE", "FG_PULMEMB", "FG_AORTANEUR",
         "FG_AF", "FG_VHD", "FG_NONISCH_CM"],
        ["varicose veins", "venous thromboembolism", "pulmonary embolism", "aortic aneurysm",
         "atrial fibrillation", "valvular heart disease", "non-ischaemic cardiomyopathy"],
        panel)
    forest(ax, items, rules=[(2.5, "cardiac anchors")])
    ax.set_xlabel("effect of eBMD (log OR per SD)", fontsize=6.6)
    ax.set_title("Cardiac, not systemic matrix", loc="left", fontsize=7.6)


def panel_b(ax, panel):
    items = rows_for(
        ["FG_OTHARR", "FG_PAROXTAC", "FG_CONDUCTIO", "FG_AVBLOCK", "FG_AF_REIMB", "FG_AF"],
        ["other arrhythmias", "paroxysmal tachycardia", "conduction disorders", "AV block",
         "atrial fibrillation (reimbursement)", "atrial fibrillation"],
        panel)
    forest(ax, items, rules=[(1.5, "atrial")])
    ax.set_xlabel("effect of eBMD (log OR per SD)", fontsize=6.6)
    ax.set_title("Atrial fibrillation, not the conduction system", loc="left", fontsize=7.6)


def panel_c(ax, panel):
    """The three aggregates against the two specific endpoints, with the pooled contrast.

    Strict monotonicity fails on the three aggregates, whose estimates all sit within one
    standard error of zero and differ from each other in sign. The test that carries
    information is the pooled contrast, so that is what the panel draws.
    """
    items = rows_for(
        ["FG_CVD_ALL", "FG_IHD", "FG_CHD", "FG_HF", "FG_NONISCH_CM"],
        ["cardiovascular disease, all", "ischaemic heart disease", "major CHD event",
         "heart failure", "non-ischaemic cardiomyopathy"],
        panel)
    forest(ax, items, rules=[(1.5, "specific")])
    ax.set_xlabel("effect of eBMD (log OR per SD)", fontsize=6.6)
    ax.set_title("Aggregates carry no effect", loc="left", fontsize=7.6)

    agg = [s for s in items if s["key"] in ("FG_CVD_ALL", "FG_IHD", "FG_CHD")]
    spec = [s for s in items if s["key"] in ("FG_HF", "FG_NONISCH_CM")]
    if len(agg) >= 2 and len(spec) >= 2:
        def pool(group):
            b = np.array([s["b"] for s in group])
            se = np.array([(s["hi"] - s["lo"]) / 3.92 for s in group])
            w = 1 / se ** 2
            return float(np.sum(w * b) / np.sum(w)), float(np.sqrt(1 / np.sum(w)))
        ba, sea = pool(agg)
        bs, ses = pool(spec)
        diff = bs - ba
        sed = float(np.sqrt(sea ** 2 + ses ** 2))
        p_val = float(2 * stats.norm.sf(abs(diff / sed)))
        ax.annotate(f"pooled: aggregates {ba:+.3f}, specific {bs:+.3f}\n"
                    f"difference {diff:+.3f} (p = {p_val:.1e})",
                    xy=(0.02, 0.03), xycoords="axes fraction", ha="left", va="bottom",
                    fontsize=5.5, color="0.25")


def panel_d(ax):
    lad = pd.read_csv(TAB / "ladder.csv").set_index("outcome")
    avc = pd.read_csv(TAB / "avc.csv").set_index("outcome")
    items = []
    for key, label in [("FG_VHD", "valvular disease, ICD aggregate"),
                       ("FG_CAVS", "calcific aortic stenosis, operated"),
                       ("AVC", "aortic valve calcification")]:
        src = avc if key in avc.index else lad
        if key not in src.index:
            continue
        r = src.loc[key]
        items.append(dict(key=key, label=label, b=float(r["b"]), lo=float(r["lo"]),
                          hi=float(r["hi"]), p=float(r["p"]),
                          power=float(r["power_ref"]), cases=int(r["cases"])))
    forest(ax, items)
    ax.set_xlabel("effect of eBMD (log OR per SD)", fontsize=6.6)
    ax.set_title("Whether the valve effect is calcific stays open", loc="left", fontsize=7.6)


def main():
    panel = pd.read_csv(TAB / "panel.csv")

    apply_figure_style(sizes=(7.6, 6.6, 5.8))
    fig = plt.figure(figsize=(7.4, 5.0))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.06, 1.0], hspace=0.62, wspace=0.86)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    panel_a(ax_a, panel)
    panel_b(ax_b, panel)
    panel_c(ax_c, panel)
    panel_d(ax_d)

    for ax, letter in zip([ax_a, ax_b, ax_c, ax_d], "abcd"):
        panel_letter(ax, letter)

    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], color=PROTECT, marker="o", ms=4.6, lw=1.5, label="protective, p < 0.05"),
        Line2D([], [], color=NULL, marker="o", ms=4.6, lw=1.5, label="null"),
        Line2D([], [], color=OPEN, marker="o", ms=4.6, lw=1.2, mfc="white", mew=1.2,
               label="null, power < 65 %"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.055), ncol=3,
               frameon=False, fontsize=6.0, handlelength=1.7, handletextpad=0.55,
               columnspacing=1.6)

    FIG.mkdir(exist_ok=True)
    out = FIG / "fig11_discriminating_panel.png"
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
