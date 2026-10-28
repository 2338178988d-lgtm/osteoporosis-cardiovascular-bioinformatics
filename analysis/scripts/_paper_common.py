"""Shared drawing helpers for the six paper figures (figures/paper/).

Keeps the six figures on one visual grammar: one colour per meaning, one type scale, forests
drawn the same way everywhere, and a collision scan run on every render.

Colour bindings (do not reassign):
  PROTECT  blue    a protective (negative) estimate that is significant
  NULL     grey    an estimate compatible with no effect
  WEAK     orange  a null that could not have detected the effect of interest
  POOL     black   a pooled / anchor estimate
  FI       blue    Finnish cohorts in the replication panels (same blue as PROTECT: the Finnish
                   cohorts are the protective ones, so one hue carries one idea)
  OTHER    grey    non-Finnish cohorts
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = None  # set by each figure script before use

PROTECT = "#1f6fb4"
NULL = "#8a8a8a"
WEAK = "#b45f1f"
POOL = "#111111"
FI = PROTECT
OTHER = NULL


def load_style():
    """Apply the project figure style. Call once, before the first plotting call.

    Loads the vendored copy first (`analysis/scripts/_figure_style_runtime.py`, identical to
    the skill runtime that generated it), so the figures rebuild on any machine and do not
    depend on an editor-specific skill directory being present.
    """
    here = Path(__file__).resolve().parent
    candidate = here / "_figure_style_runtime.py"
    if candidate.exists():
        exec(compile(open(candidate, encoding="utf-8").read(), "runtime.py", "exec"),
             globals())
        return candidate
    raise FileNotFoundError(
        "figure-style runtime not found; expected analysis/scripts/_figure_style_runtime.py")


def forest(ax, items, xlabel="", zero_line=True, show_power=False, label_size=6.0):
    """One row per estimate.

    items: list of dicts with label, b, lo, hi, and optional p, power, kind
    kind: 'protect' | 'null' | 'weak' | 'pool'
    """
    y = np.arange(len(items))[::-1]
    if zero_line:
        ax.axvline(0, color="0.45", lw=0.9, zorder=1)
    for yy, s in zip(y, items):
        kind = s.get("kind", "null")
        colour = {"protect": PROTECT, "null": NULL, "weak": WEAK, "pool": POOL}[kind]
        hollow = kind == "weak"
        ax.plot([s["lo"], s["hi"]], [yy, yy], color=colour,
                lw=2.0 if kind == "pool" else 1.5,
                solid_capstyle="butt", zorder=2, alpha=0.6 if hollow else 1.0)
        marker = "D" if kind == "pool" else "o"
        ax.plot(s["b"], yy, marker, color="white" if hollow else colour,
                ms=5.8 if kind == "pool" else 4.8,
                mec=colour, mew=1.2 if hollow else 0.55, zorder=3)
        ax.annotate(f"{s['b']:+.3f}", (s["hi"], yy), textcoords="offset points",
                    xytext=(4, -1), va="center", fontsize=5.9, color=colour)
        if show_power and "power" in s and np.isfinite(s["power"]):
            ax.annotate(f"{s['power']:.0%}", (1.0, yy),
                        xycoords=("axes fraction", "data"), textcoords="offset points",
                        xytext=(26, -1), va="center", ha="right", fontsize=5.3,
                        color="0.5", annotation_clip=False)
    ax.set_yticks(y)
    ax.set_yticklabels([s["label"] for s in items], fontsize=label_size)
    ax.set_ylim(-0.9, len(items) - 0.4)
    ax.margins(x=0.22)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=6.7)
    if show_power:
        ax.annotate("power", (1.0, len(items) - 0.75), xycoords=("axes fraction", "data"),
                    textcoords="offset points", xytext=(26, 0), ha="right",
                    fontsize=5.3, color="0.5", annotation_clip=False)
    return y


def kind_of(b, p, power=None, thr=0.05):
    """Classify an estimate for drawing: protected, informative null, or underpowered null."""
    if p < thr and b < 0:
        return "protect"
    if power is not None and np.isfinite(power) and power < 0.65:
        return "weak"
    return "null"


def block_label(ax, y0, y1, text):
    """Italic group header spanning rows y0..y1 (data coordinates)."""
    ax.annotate(text, xy=(0.0, y1), xycoords=("axes fraction", "data"),
                textcoords="offset points", xytext=(-2, 7), ha="left", va="bottom",
                fontsize=5.8, color="0.4", style="italic")


def collide(fig):
    """Assert no visible text overlaps another or a spine, and all text is inside the canvas."""
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


def save(fig, name, root):
    """Save to figures/paper/ and report the collision scan."""
    out = root / "figures" / "paper" / name
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    hits = collide(fig)
    print(f"saved {out.relative_to(root).as_posix()}")
    print(f"collision scan: {len(hits)} overlaps")
    for h in hits[:8]:
        print("   ", ascii(h))
    return hits

