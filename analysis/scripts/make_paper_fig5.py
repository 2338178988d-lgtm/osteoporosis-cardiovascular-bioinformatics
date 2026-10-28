"""Figure 5 - in the population, low bone mass predicts death but is not cross-sectionally
associated with cardiovascular disease.

Panel a  prevalent cardiovascular disease by bone density: the crude association disappears once
         age, sex, lifestyle and metabolic factors are in the model
Panel b  mortality by bone density: fully adjusted, it does not disappear
Panel c  the four cardiovascular outcomes compared against each other: none differs from the
         others

Claim: "Low bone mass is a prognostic marker, not a comorbid state - and the population layer
never had an outcome-specific signal to begin with." Panel c is the part that changes with this
project: it was written as a limitation, and once the genetic layer showed the effect spanning
several non-coronary endpoints rather than only heart failure, it became support.

All three panels are vertical forests so the four outcomes can be read across them.

Run from the project root:  python analysis/scripts/make_paper_fig5.py
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
from matplotlib.lines import Line2D

TAB = ROOT / "results" / "tables" / "01-nhanes"
OUTCOMES = [("chf", "heart failure"), ("chd", "coronary heart disease"),
            ("mi", "myocardial infarction"), ("stroke", "stroke")]
SHORT = {"chf": "heart failure", "chd": "coronary disease",
         "mi": "myocardial infarction", "stroke": "stroke"}


def panel_a(ax):
    """Prevalent disease: crude versus adjusted, one row per outcome."""
    d = pd.read_csv(TAB / "table12_outcome_comparison.csv").set_index("model")
    y = np.arange(len(OUTCOMES))[::-1]
    ax.axvline(1.0, color="0.5", lw=0.8, zorder=1)

    for yy, (out, _) in zip(y, OUTCOMES):
        for off, model, colour in [(0.17, "M1", pc.PROTECT), (-0.17, "M4", pc.POOL)]:
            or_ = float(d.loc[model, f"or.{out}"])
            lo, hi = float(d.loc[model, f"lo.{out}"]), float(d.loc[model, f"hi.{out}"])
            ax.plot([lo, hi], [yy + off, yy + off], color=colour, lw=1.5,
                    solid_capstyle="butt", zorder=2)
            ax.plot(or_, yy + off, "o", color=colour, ms=4.6, mec="white", mew=0.5, zorder=3)
        ax.annotate(f"{float(d.loc['M1', f'or.{out}']):.2f}",
                    (float(d.loc["M1", f"hi.{out}"]), yy + 0.17),
                    textcoords="offset points", xytext=(3, 0), ha="left", va="center",
                    fontsize=5.4, color=pc.PROTECT)
        ax.annotate(f"{float(d.loc['M4', f'or.{out}']):.2f}",
                    (float(d.loc["M4", f"hi.{out}"]), yy - 0.17),
                    textcoords="offset points", xytext=(3, 0), ha="left", va="center",
                    fontsize=5.4, color=pc.POOL)

    ax.set_yticks(y)
    ax.set_yticklabels([SHORT[o] for o, _ in OUTCOMES], fontsize=6.0)
    ax.set_xlim(0.62, 1.80)
    ax.set_ylim(-0.75, len(OUTCOMES) - 0.25)
    ax.set_xlabel("odds ratio, osteoporosis vs normal BMD", fontsize=6.6)
    ax.set_xticks([0.7, 1.0, 1.3, 1.6])
    ax.set_title("Prevalent disease: three of four null", loc="left", fontsize=7.4)
    handles = [Line2D([], [], color=pc.PROTECT, marker="o", ms=4.4, lw=1.5, label="crude"),
               Line2D([], [], color=pc.POOL, marker="o", ms=4.4, lw=1.5,
                      label="fully adjusted")]
    ax.legend(handles=handles, frameon=False, fontsize=5.5, loc="lower right",
              handlelength=1.4, handletextpad=0.5, borderpad=0.1, labelspacing=0.3)


def panel_b(ax):
    """Mortality by bone-density category, fully adjusted."""
    d = pd.read_csv(TAB / "table4_survey_cox_mortality.csv")
    d = d[(d["model"].str.startswith("M3")) & (d["subgroup"] == "All")]
    rows = [("dead", "Osteoporosis", "all-cause death"),
            ("dead", "Osteopenia", "all-cause death"),
            ("cvd_death_num", "Osteoporosis", "cardiovascular death"),
            ("cvd_death_num", "Osteopenia", "cardiovascular death")]

    y = np.arange(len(rows))[::-1]
    ax.axvline(1.0, color="0.5", lw=0.8, zorder=1)
    for yy, (outcome, cat, _) in zip(y, rows):
        r = d[(d["outcome"] == outcome) & (d["term"] == f"bmd_cat{cat}")].iloc[0]
        colour = pc.PROTECT if r["p"] < 0.05 else pc.NULL
        ax.plot([r["lcl"], r["ucl"]], [yy, yy], color=colour, lw=1.6,
                solid_capstyle="butt", zorder=2)
        ax.plot(r["HR"], yy, "o", color=colour, ms=5.0, mec="white", mew=0.5, zorder=3)
        ax.annotate(f"{r['HR']:.2f}", (r["ucl"], yy), textcoords="offset points",
                    xytext=(4, 0), va="center", fontsize=5.6, color=colour)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{cat.lower()}, {name.replace('cardiovascular', 'CV')}" for _, cat, name in rows], fontsize=5.7)
    ax.set_xlim(0.7, 2.9)
    ax.set_ylim(-0.75, len(rows) - 0.25)
    ax.set_xlabel("hazard ratio (fully adjusted)", fontsize=6.6)
    ax.set_xticks([1.0, 1.5, 2.0, 2.5])
    ax.set_title("Incident death: association holds", loc="left", fontsize=7.4)


def panel_c(ax):
    """The four outcomes are not different from each other, per SD of bone density."""
    d = pd.read_csv(TAB / "table12_outcome_comparison.csv").set_index("model")
    het = pd.read_csv(TAB / "table13_outcome_heterogeneity.csv")

    y = np.arange(len(OUTCOMES))[::-1]
    ax.axvline(1.0, color="0.5", lw=0.8, zorder=1)
    for yy, (out, _) in zip(y, OUTCOMES):
        or_ = float(d.loc["M4", f"or.{out}"])
        lo, hi = float(d.loc["M4", f"lo.{out}"]), float(d.loc["M4", f"hi.{out}"])
        ax.plot([lo, hi], [yy, yy], color=pc.NULL, lw=1.6,
                solid_capstyle="butt", zorder=2)
        ax.plot(or_, yy, "o", color=pc.NULL, ms=5.0, mec="white", mew=0.5, zorder=3)
        ax.annotate(f"{or_:.2f}", (hi, yy), textcoords="offset points",
                    xytext=(4, 0), va="center", fontsize=5.6, color=pc.NULL)
    ax.set_yticks(y)
    ax.set_yticklabels([SHORT[o] for o, _ in OUTCOMES], fontsize=6.0)
    ax.set_xlim(0.80, 1.40)
    ax.set_xticks([0.9, 1.0, 1.1, 1.2, 1.3])
    ax.set_ylim(-0.75, len(OUTCOMES) - 0.25)
    ax.set_xlabel("odds ratio per SD lower bone density (adjusted)", fontsize=6.6)
    ax.set_title("No outcome differs from the others", loc="left", fontsize=7.4)
    pairs = ", ".join(f"{r['comparison'].replace('CHF', 'HF').replace('chd', 'CHD').replace('mi', 'MI')}"
                      f" {r['p']:.2f}" for _, r in het.iterrows())
    ax.annotate(f"heterogeneity vs heart failure: {pairs}\nstroke remains associated; "
                f"it does not differ from heart failure", xy=(0.5, -0.36),
                xycoords="axes fraction", ha="center", va="top", fontsize=5.0,
                color="0.35")


def main():
    pc.apply_figure_style(sizes=(7.6, 6.6, 5.8))
    fig = plt.figure(figsize=(7.6, 2.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.05, 1.0], wspace=1.05)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(3))

    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)

    for ax, letter in zip([ax_a, ax_b, ax_c], "abc"):
        pc.panel_letter(ax, letter)

    pc.save(fig, "fig5_population.png", ROOT)
    return 0


if __name__ == "__main__":
    sys.exit(main())

