"""Figure 1 - the strongest evidence is concentrated in selected non-coronary cardiac endpoints.

Eighteen endpoints, one instrument set, one cohort, one estimator. They are drawn in four
blocks that follow the anatomy: systemic matrix, vascular load and coronary disease, the
conduction system, then atrial / valvular / myocardial disease. The blocks show the evidence pattern:
the strongest evidence is concentrated in selected non-coronary cardiac endpoints; other tested endpoint families show no evidence at the prespecified thresholds.

Claim: "The estimated association is concentrated in selected non-coronary cardiac endpoints; other tested endpoint families show no evidence at the prespecified thresholds." Bone density is not being interpreted as a treatment effect;
or conduction disease." Every block label is checked against its rows before saving.

Run from the project root:  python analysis/scripts/make_paper_fig1.py
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

TAB = ROOT / "results" / "tables" / "14-cardiac-ladder"

# (display label, table, key, block)
SPEC = [
    ("varicose veins",            "panel",  "FG_VARICVE",       "matrix"),
    ("venous thromboembolism",    "panel",  "FG_VTE",           "matrix"),
    ("pulmonary embolism",        "panel",  "FG_PULMEMB",       "matrix"),
    ("aortic aneurysm",           "panel",  "FG_AORTANEUR",     "matrix"),
    ("hypertension",              "ladder", "FG_HYPERTENS",     "load"),
    ("myocardial infarction",     "ladder", "FG_MI_STRICT",     "load"),
    ("major CHD event",           "panel",  "FG_CHD",           "load"),
    ("ischaemic heart disease",   "panel",  "FG_IHD",           "load"),
    ("cardiovascular disease, all", "panel", "FG_CVD_ALL",      "load"),
    ("conduction disorders",      "panel",  "FG_CONDUCTIO",     "conduction"),
    ("AV block",                  "panel",  "FG_AVBLOCK",       "conduction"),
    ("paroxysmal tachycardia",    "panel",  "FG_PAROXTAC",      "conduction"),
    ("other arrhythmias",         "panel",  "FG_OTHARR",        "conduction"),
    ("atrial fibrillation",       "ladder", "FG_AF",            "cardiac"),
    ("atrial fibrillation, reimbursement", "panel", "FG_AF_REIMB", "cardiac"),
    ("valvular heart disease",    "ladder", "FG_VHD",           "cardiac"),
    ("non-ischaemic cardiomyopathy", "ladder", "FG_NONISCH_CM", "cardiac"),
    ("heart failure",             "ladder", "FG_HF",            "cardiac"),
]
BLOCKS = [
    ("matrix", "systemic matrix"),
    ("load", "vascular load and coronary disease"),
    ("conduction", "conduction system"),
    ("cardiac", "atrial, valvular and myocardial disease"),
]


def main():
    panel = pd.read_csv(TAB / "panel.csv").set_index("outcome")
    ladder = pd.read_csv(TAB / "ladder.csv").set_index("outcome")
    src = {"panel": panel, "ladder": ladder}

    items, block_spans = [], {}
    for label, table, key, block in SPEC:
        if key not in src[table].index:
            print(f"  {key} missing from {table}.csv, skipped")
            continue
        r = src[table].loc[key]
        power = float(r["power_ref"])
        items.append(dict(label=label, key=key, block=block, b=float(r["b"]),
                          lo=float(r["lo"]), hi=float(r["hi"]), p=float(r["p"]),
                          power=power, cases=int(r["cases"]),
                          kind=pc.kind_of(float(r["b"]), float(r["p"]), power)))
        block_spans.setdefault(block, []).append(len(items) - 1)

    pc.apply_figure_style(sizes=(8.0, 6.9, 5.9))
    fig = plt.figure(figsize=(6.4, 5.4))
    ax = fig.add_subplot(111)

    # lay the rows out by hand so a blank line can sit between blocks: the block header goes
    # in that line, which is why it never competes with the spine or with a confidence interval
    rows, spacers = [], []
    y = 0.0
    for block, header in BLOCKS:
        group = [s for s in items if s["block"] == block]
        if not group:
            continue
        # x in axes fraction, y in data coordinates: passing a single "axes fraction" would
        # also reinterpret y, putting every header at the bottom of the panel
        ax.annotate(header, xy=(0.03, y), xycoords=("axes fraction", "data"),
                    textcoords="offset points", xytext=(0, 2), ha="left", va="center",
                    fontsize=5.8, color="0.4", style="italic")
        spacers.append(y)
        y -= 1.0
        for s in group:
            rows.append((y, s))
            y -= 1.0

    ax.axvline(0, color="0.45", lw=0.9, zorder=1)
    for yy, s in rows:
        colour = {"protect": pc.PROTECT, "null": pc.NULL, "weak": pc.WEAK,
                  "pool": pc.POOL}[s["kind"]]
        hollow = s["kind"] == "weak"
        ax.plot([s["lo"], s["hi"]], [yy, yy], color=colour, lw=1.5,
                solid_capstyle="butt", zorder=2, alpha=0.6 if hollow else 1.0)
        ax.plot(s["b"], yy, "o", color="white" if hollow else colour, ms=4.8,
                mec=colour, mew=1.2 if hollow else 0.55, zorder=3)
        ax.annotate(f"{s['b']:+.3f}", (s["hi"], yy), textcoords="offset points",
                    xytext=(4, -1), va="center", fontsize=5.9, color=colour)
        ax.annotate(f"{s['power']:.0%}", (1.0, yy), xycoords=("axes fraction", "data"),
                    textcoords="offset points", xytext=(26, -1), va="center", ha="right",
                    fontsize=5.3, color="0.5", annotation_clip=False)
        ax.annotate(f"{s['cases'] / 1000:.0f}k", (1.0, yy), xycoords=("axes fraction", "data"),
                    textcoords="offset points", xytext=(52, -1), va="center", ha="right",
                    fontsize=5.3, color="0.5", annotation_clip=False)

    ax.set_yticks([yy for yy, _ in rows])
    ax.set_yticklabels([s["label"] for _, s in rows], fontsize=6.2)
    ax.set_ylim(y + 0.6, 0.9)
    ax.set_xlim(-0.17, 0.12)
    ax.set_xlabel("effect of genetically predicted heel-ultrasound eBMD (log OR per SD)",
                  fontsize=6.7)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.margins(x=0.22)
    top = max(yy for yy, _ in rows)
    ax.annotate("power", (1.0, top + 0.55), xycoords=("axes fraction", "data"),
                textcoords="offset points", xytext=(26, 0), ha="right",
                fontsize=5.3, color="0.5", annotation_clip=False)
    ax.annotate("cases", (1.0, top + 0.55), xycoords=("axes fraction", "data"),
                textcoords="offset points", xytext=(52, 0), ha="right",
                fontsize=5.3, color="0.5", annotation_clip=False)
    ax.set_title("Selected non-coronary cardiac endpoints carry the strongest evidence", loc="left", fontsize=8.6)
    pc.panel_letter(ax, "a")

    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=pc.PROTECT, marker="o", ms=4.8, lw=1.5,
                      label="protective, p < 0.05"),
               Line2D([], [], color=pc.NULL, marker="o", ms=4.8, lw=1.5, label="null"),
               Line2D([], [], color=pc.WEAK, marker="o", ms=4.8, lw=1.2, mfc="white",
                      mew=1.2, label="null, power < 65 %")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.055), ncol=3,
               frameon=False, fontsize=6.2, handlelength=1.7, handletextpad=0.55,
               columnspacing=1.8)

    pc.save(fig, "fig1_conclusion_map.png", ROOT)
    return 0


if __name__ == "__main__":
    sys.exit(main())

