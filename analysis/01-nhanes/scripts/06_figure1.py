"""Figure 1 - NHANES bone-vascular association: four panels for one claim
("lower femoral neck BMD tracks a higher burden of cardiovascular disease and
of subsequent death in US adults aged 40+").

Panels
  a  survey-weighted prevalence of CVD by WHO BMD category, both sexes and by sex
  b  adjusted ORs for prevalent CVD (forest, four nested models)
  c  dose-response: adjusted OR across the T-score continuum (RCS)
  d  weighted Kaplan-Meier survival by BMD category

Inputs
  results/tables/01-nhanes/table2_cvd_prevalence_by_bmd.csv
  results/tables/01-nhanes/table3_survey_logistic_cvd.csv
  results/tables/01-nhanes/table8_rcs_cvd_by_t_score.csv
  results/tables/01-nhanes/table9_km_allcause_by_bmd.csv
Output
  figures/fig1_nhanes_bmd_cvd.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from _figure_style_runtime import apply_figure_style

apply_figure_style(frame="open", sizes=(8, 7, 6))
T = "results/tables/01-nhanes/"
CAT = ["Normal", "Osteopenia", "Osteoporosis"]
# one hue per BMD category, held constant across all panels
COL = {"Normal": "#4C72B0", "Osteopenia": "#DD8452", "Osteoporosis": "#B03A5B"}

fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4))
(axa, axb), (axc, axd) = axes

# ---------------------------------------------------------------- panel a
prev = pd.read_csv(T + "table2_cvd_prevalence_by_bmd.csv")
prev["bmd_cat"] = pd.Categorical(prev["bmd_cat"], CAT, ordered=True)
both = prev[prev.sex == "Both"].sort_values("bmd_cat")
xs = np.arange(len(CAT))
bars = axa.bar(xs, both.weighted_prevalence * 100, width=0.58,
               color=[COL[c] for c in CAT], edgecolor="black", linewidth=0.6)
for x, (_, r) in zip(xs, both.iterrows()):
    axa.errorbar(x, r.weighted_prevalence * 100,
                 yerr=[[r.weighted_prevalence * 100 - r.ci_l * 100],
                       [r.ci_u * 100 - r.weighted_prevalence * 100]],
                 fmt="none", ecolor="black", elinewidth=0.9, capsize=2.5)
    axa.text(x, r.ci_u * 100 + 0.9, f"{r.weighted_prevalence*100:.1f}%",
             ha="center", va="bottom", fontsize=7)
# sex-specific points, offset, open markers so they never read as extra bars
for sex, dx, mk in (("Female", -0.16, "o"), ("Male", 0.16, "s")):
    s = prev[prev.sex == sex].set_index("bmd_cat").loc[CAT]
    axa.plot(xs + dx, s.weighted_prevalence * 100, mk, ms=4.2, mfc="white",
             mec="black", mew=0.8, ls="none", zorder=5)
axa.set_xticks(xs); axa.set_xticklabels(CAT)
axa.set_ylabel("Prevalent CVD (%)")
axa.set_ylim(0, 32)
axa.set_title("CVD prevalence rises across WHO BMD strata", loc="left")
axa.text(0.02, 0.96, "filled = both sexes, 95% CI\nopen 鈼?women  鈻?men",
         transform=axa.transAxes, fontsize=6, va="top")
axa.margins(0.04)

# ---------------------------------------------------------------- panel b
lg = pd.read_csv(T + "table3_survey_logistic_cvd.csv")
sub = lg[(lg.outcome == "cvd") & (lg.exposure == "bmd_cat") & (lg.subgroup == "Both")].copy()
if sub.empty:
    sub = lg[(lg.outcome == "cvd") & (lg.exposure == "bmd_cat") & (lg.subgroup == "All")].copy()
sub["grp"] = sub["term"].str.replace("bmd_cat", "", regex=False)
models = [m for m in ["M1 crude", "M2 +age/sex/race/cycle", "M3 +lifestyle/metabolic",
                      "M4 +mineral metabolism"] if m in set(sub.model)]
models = models[::-1]
ys, labels, colors = [], [], []
row = 0
for m in models:
    part = sub[sub.model == m].set_index("grp").reindex(["Osteopenia", "Osteoporosis"])
    for g in ["Osteopenia", "Osteoporosis"]:
        if g not in part.index or pd.isna(part.loc[g, "OR"]):
            continue
        r = part.loc[g]
        ys.append(row); labels.append((m, g)); colors.append(COL[g])
        axb.plot([r.lcl, r.ucl], [row, row], color=COL[g], lw=0.9, solid_capstyle="butt")
        axb.plot(r.OR, row, "o", ms=3.6, color=COL[g], mec="white", mew=0.4)
        row += 1
    row += 0.6
axb.axvline(1, color="0.35", lw=0.8, ls=(0, (2, 2)))
axb.set_xscale("log")
# explicit ticks inside the view range: leaving matplotlib's default log ticks in
# place draws labels for 0.1 and 10 that sit outside the limits
axb.set_xticks([0.6, 0.8, 1.0, 1.5, 2.0, 3.0])
axb.set_xticklabels(["0.6", "0.8", "1.0", "1.5", "2.0", "3.0"])
axb.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
axb.set_yticks([y for y, (m, g) in zip(ys, labels)])
axb.set_yticklabels([f"{g[:5].lower()}." for (m, g) in labels], fontsize=6)
# one left-hand model label per block
seen = {}
for y, (m, g) in zip(ys, labels):
    seen.setdefault(m, []).append(y)
for m, yy in seen.items():
    axb.text(-0.02, np.mean(yy), m, transform=axb.get_yaxis_transform(),
             ha="right", va="center", fontsize=6, color="0.25")
axb.set_ylabel("")
axb.set_xlabel("Adjusted odds ratio for prevalent CVD (log scale)")
axb.set_title("Association survives adjustment", loc="left")
axb.set_xlim(0.55, 3.2)
axb.margins(x=0.02)
axb.grid(axis="y", visible=False)

# ---------------------------------------------------------------- panel c
rcs = pd.read_csv(T + "table8_rcs_cvd_by_t_score.csv")
axc.fill_between(rcs.t_score, rcs.lcl, rcs.ucl, color="#4C72B0", alpha=0.18, lw=0)
axc.plot(rcs.t_score, rcs.OR, color="#4C72B0", lw=1.4)
axc.axhline(1, color="0.35", lw=0.8, ls=(0, (2, 2)))
for t, lab in [(-2.5, "osteoporosis\nthreshold"), (-1.0, "osteopenia\nthreshold")]:
    axc.axvline(t, color="0.6", lw=0.6, ls=(0, (1, 3)))
    axc.text(t - 0.08, axc.get_ylim()[1] * 0.97, lab, rotation=90, ha="right", va="top",
             fontsize=5.5, color="0.35")
axc.set_xlabel("Femoral neck T-score (sex-specific)")
axc.set_ylabel("Adjusted OR for prevalent CVD")
axc.set_title("Smooth gradient, no threshold", loc="left")
axc.set_xlim(-4.5, 2.5)
axc.margins(0.04)

# ---------------------------------------------------------------- panel d
km = pd.read_csv(T + "table9_km_allcause_by_bmd.csv")
for g in CAT:
    k = km[km.stratum == g]
    axd.plot(k.time, k.surv * 100, color=COL[g], lw=1.5)
    axd.fill_between(k.time, k.lcl * 100, k.ucl * 100, color=COL[g], alpha=0.15, lw=0)
    j = k[k.time == 14]
    if len(j):
        axd.text(14.15, float(j.surv.iloc[0]) * 100, g, color=COL[g], fontsize=6.5,
                 va="center", ha="left")
axd.set_xlim(0, 17.4)
axd.set_ylim(30, 102)
axd.set_xlabel("Follow-up from examination (years)")
axd.set_ylabel("Survival (%)")
axd.set_title("Lower BMD, steeper mortality", loc="left")
axd.margins(0.04)

import matplotlib.ticker as mt
axd.yaxis.set_major_locator(mt.MultipleLocator(20))

for ax, letter in zip([axa, axb, axc, axd], "abcd"):
    panel_letter(ax, letter)

# give axis labels clear air from their tick labels (the style ladder uses small
# fonts, and matplotlib's default labelpad leaves the boxes touching)
for ax in fig.axes:
    ax.xaxis.labelpad = 6
    ax.yaxis.labelpad = 7
    ax.tick_params(axis="both", pad=2.5)

fig.tight_layout(w_pad=1.6, h_pad=1.4)
out = "figures/fig1_nhanes_bmd_cvd.png"
import os
os.makedirs("figures", exist_ok=True)
# Scan on the laid-out canvas *before* the tight-bbox save: savefig re-derives the
# figure box at 300 dpi, so post-save bounding boxes no longer match the layout and
# produce false collisions.
fig.canvas.draw()
import matplotlib as mpl
rend = fig.canvas.get_renderer()


def _inside_view(t):
    """matplotlib keeps Text objects for axis ticks that fall outside the current
    view limits; they are clipped when drawn, so they must not count as overlaps."""
    ax = t.axes
    if ax is None:
        return True
    if t in ax.get_xticklabels(which="both"):
        lo, hi = sorted(ax.get_xlim())
        pos = t.get_position()[0]
        return lo - 1e-9 <= pos <= hi + 1e-9
    if t in ax.get_yticklabels(which="both"):
        lo, hi = sorted(ax.get_ylim())
        pos = t.get_position()[1]
        return lo - 1e-9 <= pos <= hi + 1e-9
    return True


labels_ = [(t, t.get_window_extent(rend)) for t in fig.findobj(mpl.text.Text)
           if t.get_text().strip() and t.get_visible() and _inside_view(t)]
frames = [(s, s.get_window_extent(rend)) for ax in fig.axes
          for s in ax.spines.values() if s.get_visible()]
own_ticks = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
             for ax in fig.axes}
hits = [(a.get_text(), b.get_text()) for i, (a, ba) in enumerate(labels_)
        for b, bb in labels_[i + 1:] if ba.overlaps(bb)]
hits += [(t.get_text(), "spine") for t, bt in labels_ for s, bs in frames
         if bt.overlaps(bs) and t not in own_ticks[s.axes]]
print("collisions:", hits)
assert not hits, hits

fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
print("saved", out)


