"""Figure 3 - the effect is diffuse at every level, and no annotated set carries it.

Panel a  posterior probability of a shared causal variant, one value per instrument locus out of
         600; exactly one locus reaches PP.H4 > 0.8 and removing it does not move the estimate
Panel b  the seven pre-specified bone pathways: the "unassigned" instruments carry the whole
         effect and the classic bone genes carry the opposite sign
Panel c  the genome-wide gene-set scan: the largest observed |z| sits inside the null
         distribution, in both cohorts

Claim: "No single locus, pathway or annotated gene set carries this effect." Each panel is one
level of that claim; none of them is a supplementary check.

Run from the project root:  python analysis/scripts/make_paper_fig3.py
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

TAB = ROOT / "results" / "tables"

GROUP_LABEL = {
    "unassigned": "unassigned (561 instruments)",
    "wnt_osteoblast": "WNT / osteoblast (19)",
    "osteoclast_resorption": "osteoclast resorption (4)",
    "sex_hormone": "sex hormone (2)",
    "height_growth": "height / growth (2)",
    "immune_inflammatory": "immune / inflammatory (2)",
    "matrix_mineralisation": "matrix mineralisation (1)",
    "adiposity_bodycomp": "adiposity / body composition (1)",
}


def panel_a(ax):
    """PP.H4 across 600 loci: one locus stands out and it is not doing the work."""
    d = pd.read_csv(TAB / "09-integration" / "coloc_loci.csv")
    pp = d["PP_H4"].clip(lower=1e-4)

    bins = np.logspace(np.log10(1e-4), np.log10(1.0), 34)
    ax.hist(pp, bins=bins, color=pc.NULL, edgecolor="white", linewidth=0.4)
    ax.set_xscale("log")
    ax.set_xlabel("posterior probability of a shared causal variant (PP.H4)", fontsize=6.7)
    ax.set_ylabel("instrument loci", fontsize=6.7)
    ax.set_yscale("log")
    ax.set_ylim(0.8, 1400)

    above = d[d["PP_H4"] >= 0.8]
    ax.annotate(f"1 of {len(d)} loci above 0.8", xy=(0.04, 0.96), xycoords="axes fraction",
                ha="left", va="top", fontsize=5.9, color=pc.PROTECT)
    ax.annotate("median PP.H4 = 0.002", xy=(0.04, 0.855), xycoords="axes fraction",
                ha="left", va="top", fontsize=5.7, color="0.35")
    ax.annotate("removing that locus:\n-0.0640 to -0.0639", xy=(0.04, 0.75),
                xycoords="axes fraction", ha="left", va="top", fontsize=5.5, color="0.35")
    ax.set_title("One shared locus in 600", loc="left", fontsize=8.0)
    if len(above):
        ax.plot([above["PP_H4"].iloc[0]] * 2, [0.8, 3], color=pc.PROTECT, lw=1.6, zorder=4)
        ax.annotate("ABO", (above["PP_H4"].iloc[0], 3.4), ha="right", va="bottom",
                    fontsize=5.6, color=pc.PROTECT,
                    xytext=(-4, 0), textcoords="offset points")
    return above


def panel_b(ax):
    """Contribution of each pre-specified pathway group to the heart-failure estimate."""
    d = pd.read_csv(TAB / "12-pathway-stratified" / "pathway_contributions.csv")
    d = d[d["outcome"] == "FG_HF"].copy()
    d = d.reindex(d["contribution_b"].abs().sort_values().index)

    y = np.arange(len(d))
    ax.axvline(0, color="0.45", lw=0.9, zorder=1)
    for yy, (_, r) in zip(y, d.iterrows()):
        unassigned = r["group"] == "unassigned"
        colour = pc.POOL if unassigned else pc.NULL
        ax.plot([0, r["contribution_b"]], [yy, yy], color=colour,
                lw=2.4 if unassigned else 1.3, solid_capstyle="butt", zorder=2)
        ax.plot(r["contribution_b"], yy, "o", color=colour, ms=5.4 if unassigned else 4.2,
                mec="white", mew=0.5, zorder=3)
        pct = r["share_pct"]
        pct_txt = f"{pct:+.0f} %" if abs(pct) >= 0.5 else "0 %"
        ax.annotate(pct_txt, (0.0, yy), textcoords="offset points",
                    xytext=(4, -1) if r["contribution_b"] < 0 else (-4, -1),
                    ha="left" if r["contribution_b"] < 0 else "right",
                    va="center", fontsize=5.7, color="0.25")
    ax.set_yticks(y)
    ax.set_yticklabels([GROUP_LABEL.get(g, g) for g in d["group"]], fontsize=6.0)
    ax.set_xlabel("contribution to the pooled estimate (log OR)", fontsize=6.7)
    ax.set_ylim(-1.6, len(d) - 0.4)
    lo = float(d["contribution_b"].min()) - 0.026
    hi = float(d["contribution_b"].max()) + 0.012
    ax.set_xlim(lo, hi)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    ax.set_title("The unassigned instruments carry all of it", loc="left", fontsize=8.0)
    ax.annotate("classic bone genes carry\nthe opposite sign", xy=(0.98, 0.04),
                xycoords="axes fraction", ha="right", va="bottom", fontsize=5.6,
                color="0.3")


def panel_c(ax):
    """Observed max |z| over 9,427 gene sets against the permutation null."""
    null = pd.read_csv(TAB / "12-pathway-stratified" / "geneset_null.csv")
    fwer = pd.read_csv(TAB / "12-pathway-stratified" / "geneset_fwer.csv").set_index("outcome")

    lo = min(null["max_abs_z"].min(), fwer["observed_max_abs_z"].min()) - 0.4
    hi = max(null["max_abs_z"].max(), fwer["observed_max_abs_z"].max()) + 0.4
    bins = np.linspace(lo, hi, 30)

    fi = fwer.loc["FG_HF"]
    he = fwer.loc["HF"]
    ax.hist(null[null["outcome"] == "FG_HF"]["max_abs_z"], bins=bins, color=pc.NULL,
            alpha=0.55, edgecolor="white", linewidth=0.3,
            label=f"FinnGen null; observed |z| = {fi['observed_max_abs_z']:.2f}, "
                  f"p = {fi['p_familywise']:.2f}")
    ax.hist(null[null["outcome"] == "HF"]["max_abs_z"], bins=bins, color="white",
            edgecolor=pc.POOL, linewidth=0.6, histtype="step",
            label=f"HERMES null; observed |z| = {he['observed_max_abs_z']:.2f}, "
                  f"p = {he['p_familywise']:.2f}")
    ax.set_xlabel("largest |z| over 9,427 gene sets", fontsize=6.7)
    ax.set_ylabel("permutations", fontsize=6.7)
    ax.set_ylim(0, None)

    for outcome, colour, style, width in [("FG_HF", "#4a4a4a", "-", 1.7),
                                          ("HF", pc.POOL, (0, (3, 2)), 1.5)]:
        r = fwer.loc[outcome]
        ax.axvline(r["observed_max_abs_z"], color=colour, lw=width, ls=style, zorder=4)
    ax.set_title("The largest gene set is not larger than chance", loc="left", fontsize=8.0)
    ax.legend(frameon=False, fontsize=5.2, loc="upper center", bbox_to_anchor=(0.5, -0.26),
              ncol=1, handlelength=1.4, handletextpad=0.5, borderpad=0.1, labelspacing=0.25)


def main():
    pc.apply_figure_style(sizes=(8.0, 6.9, 5.9))
    fig = plt.figure(figsize=(7.6, 3.1))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.30, 1.0], wspace=1.05)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))

    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)

    for ax, letter in zip([ax_a, ax_b, ax_c], "abc"):
        pc.panel_letter(ax, letter)

    pc.save(fig, "fig3_diffuse.png", ROOT)
    return 0


if __name__ == "__main__":
    sys.exit(main())

