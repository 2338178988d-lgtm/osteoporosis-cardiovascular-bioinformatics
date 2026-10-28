"""Two-tissue DEG intersection: bone side x vascular side.

The design this project has to work with is "separate discovery per tissue, then ask
whether the two sides converge", because no public dataset contains both bone and
vascular transcriptomes from the same people.

The bone-side blood datasets turn out to carry very weak differential signal (max
|log2FC| ~0.7), so an intersection defined as FDR < 0.05 AND |log2FC| > 1 is empty by
construction. This script therefore reports the overlap across a grid of thresholds and,
for each one, tests it against the overlap expected by chance (hypergeometric), so an
empty or a small overlap is a result rather than a silent zero.

Outputs (results/tables/03-bulk-deg/):
  intersection_grid.csv      overlap and chance expectation for every threshold/model
  intersection_genes.csv     direction-consistent shared genes at the primary threshold
  within_side_replication.csv reproducibility inside each tissue side

Usage:
    python analysis/03-bulk-deg/scripts/03_intersect.py
    python analysis/03-bulk-deg/scripts/03_intersect.py --primary-fdr 0.05 --primary-lfc 0.5
"""
import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "results" / "tables" / "03-bulk-deg"

BONE_DISCOVERY = ["GSE56815", "GSE56814"]
BONE_VALIDATION = ["GSE7158", "GSE7429", "GSE13850", "GSE35956"]
VASCULAR_DISCOVERY = ["GSE100927"]
VASCULAR_REPLICATION = ["GSE43292", "GSE28829"]
MECHANISM = ["GSE37558"]


def load(gse):
    f = TAB / f"deg_{gse}.csv"
    if not f.exists():
        return None
    d = pd.read_csv(f)
    return d.set_index("gene")[["logFC", "P.Value", "adj.P.Val"]]


def select(d, fdr=None, lfc=0.0, p=None):
    """Genes passing the cutoffs, split by direction."""
    m = pd.Series(True, index=d.index)
    if fdr is not None:
        m &= d["adj.P.Val"] < fdr
    if p is not None:
        m &= d["P.Value"] < p
    if lfc:
        m &= d["logFC"].abs() > lfc
    sub = d[m]
    return set(sub.index[sub["logFC"] > 0]), set(sub.index[sub["logFC"] < 0])


def overlap_test(up_a, dn_a, up_b, dn_b, universe):
    """Direction-consistent overlap and its hypergeometric expectation."""
    same_up = up_a & up_b
    same_dn = dn_a & dn_b
    observed = len(same_up) + len(same_dn)
    n_a = len(up_a) + len(dn_a)
    n_b = len(up_b) + len(dn_b)
    if universe == 0 or n_a == 0 or n_b == 0:
        return dict(observed=observed, expected=np.nan, fold=np.nan, p_hyper=np.nan)
    expected = n_a * n_b / universe
    # P(X >= observed) for drawing n_b genes from `universe` with n_a "successes"
    p_hyper = stats.hypergeom.sf(observed - 1, universe, n_a, n_b) if observed else 1.0
    return dict(observed=observed, expected=round(expected, 2),
                fold=round(observed / expected, 3) if expected else np.nan,
                p_hyper=p_hyper)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primary-fdr", type=float, default=0.05)
    ap.add_argument("--primary-lfc", type=float, default=0.5)
    a = ap.parse_args()

    tables = {g: load(g) for g in BONE_DISCOVERY + BONE_VALIDATION +
              VASCULAR_DISCOVERY + VASCULAR_REPLICATION + MECHANISM}
    tables = {k: v for k, v in tables.items() if v is not None}
    print(f"loaded {len(tables)} DEG tables: {', '.join(sorted(tables))}")

    bone_pool = set().union(*[set(t.index) for g, t in tables.items() if g in BONE_DISCOVERY + BONE_VALIDATION])
    vasc_pool = set().union(*[set(t.index) for g, t in tables.items() if g in VASCULAR_DISCOVERY + VASCULAR_REPLICATION])
    universe = len(bone_pool & vasc_pool)
    print(f"genes measurable on both sides (intersection of platforms): {universe}")

    grid = []
    models = [
        ("bone_disc_vs_vascular_disc", BONE_DISCOVERY, VASCULAR_DISCOVERY),
        ("bone_disc_vs_vascular_repl", BONE_DISCOVERY, VASCULAR_REPLICATION),
        ("bone_valid_vs_vascular_disc", BONE_VALIDATION, VASCULAR_DISCOVERY),
        ("bone_disc_vs_mechanism", BONE_DISCOVERY, MECHANISM),
        ("vascular_disc_vs_mechanism", VASCULAR_DISCOVERY, MECHANISM),
    ]
    thresholds = [dict(fdr=0.05, lfc=1.0), dict(fdr=0.05, lfc=0.5), dict(fdr=0.05, lfc=0.0),
                  dict(fdr=0.10, lfc=0.0), dict(p=0.05, lfc=0.0), dict(p=0.01, lfc=0.0)]
    for name, side_a, side_b in models:
        for thr in thresholds:
            up_a, dn_a, up_b, dn_b = set(), set(), set(), set()
            for g in side_a:
                if g in tables:
                    u, d = select(tables[g], **thr)
                    up_a |= u & (bone_pool & vasc_pool)
                    dn_a |= d & (bone_pool & vasc_pool)
            for g in side_b:
                if g in tables:
                    u, d = select(tables[g], **thr)
                    up_b |= u & (bone_pool & vasc_pool)
                    dn_b |= d & (bone_pool & vasc_pool)
            row = dict(model=name, sides="+".join(side_a) + " vs " + "+".join(side_b),
                       threshold=", ".join(f"{k}={v}" for k, v in thr.items()),
                       n_a=len(up_a) + len(dn_a), n_b=len(up_b) + len(dn_b))
            row.update(overlap_test(up_a, dn_a, up_b, dn_b, universe))
            grid.append(row)
    gdf = pd.DataFrame(grid)
    gdf.to_csv(TAB / "intersection_grid.csv", index=False, float_format="%.6g")
    print("\n" + gdf[["model", "threshold", "n_a", "n_b", "observed", "expected", "fold", "p_hyper"]]
          .to_string(index=False))

    # primary analysis: direction-consistent genes shared by both sides
    up_a, dn_a, up_b, dn_b = set(), set(), set(), set()
    for g in BONE_DISCOVERY + BONE_VALIDATION:
        if g in tables:
            u, d = select(tables[g], fdr=a.primary_fdr, lfc=a.primary_lfc)
            up_a |= u
            dn_a |= d
    for g in VASCULAR_DISCOVERY + VASCULAR_REPLICATION:
        if g in tables:
            u, d = select(tables[g], fdr=a.primary_fdr, lfc=a.primary_lfc)
            up_b |= u
            dn_b |= d
    shared = sorted((up_a & up_b) | (dn_a & dn_b))
    rows = []
    for g in shared:
        rows.append(dict(gene=g,
                         direction="up" if g in (up_a & up_b) else "down",
                         bone_logFC=round(np.mean([tables[k].loc[g, "logFC"] for k in BONE_DISCOVERY + BONE_VALIDATION
                                                   if k in tables and g in tables[k].index]), 3),
                         vascular_logFC=round(np.mean([tables[k].loc[g, "logFC"] for k in VASCULAR_DISCOVERY + VASCULAR_REPLICATION
                                                       if k in tables and g in tables[k].index]), 3)))
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "intersection_genes.csv", index=False)
    print(f"\nprimary (FDR<{a.primary_fdr}, |log2FC|>{a.primary_lfc}) direction-consistent shared genes: {len(out)}")
    if len(out):
        print(out.head(25).to_string(index=False))

    # within-side reproducibility: does either side even reproduce itself?
    rep = []
    for label, group in [("bone", BONE_DISCOVERY + BONE_VALIDATION),
                         ("vascular", VASCULAR_DISCOVERY + VASCULAR_REPLICATION)]:
        for g1, g2 in itertools.combinations([g for g in group if g in tables], 2):
            for thr in thresholds[:4]:
                u1, d1 = select(tables[g1], **thr)
                u2, d2 = select(tables[g2], **thr)
                common = len(tables[g1].index.intersection(tables[g2].index))
                r = overlap_test(u1, d1, u2, d2, common)
                rep.append(dict(side=label, pair=f"{g1}_vs_{g2}",
                                threshold=", ".join(f"{k}={v}" for k, v in thr.items()),
                                n_1=len(u1) + len(d1), n_2=len(u2) + len(d2), **r))
    rdf = pd.DataFrame(rep)
    rdf.to_csv(TAB / "within_side_replication.csv", index=False, float_format="%.6g")
    print("\nwithin-side replication (FDR<0.05, |log2FC|>0.5):")
    print(rdf[rdf["threshold"] == "fdr=0.05, lfc=0.5"]
          [["side", "pair", "n_1", "n_2", "observed", "expected", "fold", "p_hyper"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()

