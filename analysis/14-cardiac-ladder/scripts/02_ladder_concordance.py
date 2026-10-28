"""Module 14, part 2 鈥?are the protective rungs one signal or several?

The ladder found the same protective effect on atrial fibrillation, valvular heart disease and
non-ischaemic cardiomyopathy, and nothing on hypertension or myocardial infarction. Before
that is read as three confirmations, one has to ask whether the same variants carry all three
effects: the Finnish endpoints are ICD-code sets that share patients (heart failure, AF and
valve disease overlap heavily), so three concordant estimates can be one axis observed three
times. They are read as concordant, not independent.

This script answers the narrower, checkable question: at the 577 instrument variants, do the
per-variant causal estimates agree between endpoints? Per-variant Wald ratios (by / bx) are
correlated pairwise; a high correlation across all protective rungs and a low one against MI
means one axis with a distinct, non-coronary signature.

Outputs (results/tables/14-cardiac-ladder/):
  rung_concordance.csv    pairwise correlation of per-variant Wald ratios

Usage:
  python analysis/14-cardiac-ladder/scripts/02_ladder_concordance.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "analysis" / "13-robustness" / "scripts"))

from _common import mr  # noqa: E402

OUT = ROOT / "results" / "tables" / "14-cardiac-ladder"
CACHE = ROOT / "data" / "processed" / "mr" / "outcome_snps"
EXPOSURE = "eBMD"

PROTECTIVE = ["FG_HF", "FG_AF", "FG_VHD", "FG_NONISCH_CM"]
CONTRAST = ["FG_HYPERTENS", "FG_MI_STRICT"]


def ratios(inst, key):
    """Per-variant Wald ratio for one outcome, indexed by SNP.

    Read from the cache written by 01_ladder.py, not through `mr.outcome_at_snps`: that
    helper returns nothing unless the cache covers *every* requested variant, while the
    whole point here is to use the variants the cache does hold.
    """
    path = CACHE / f"{key}.csv"
    if not path.exists():
        return None
    out = pd.read_csv(path)
    if out.empty:
        return None
    d = mr.harmonise(inst, out).set_index("SNP")
    return (d["by"] / d["bx"]).rename(key)


def main():
    inst = mr.load_instruments(EXPOSURE)
    if inst is None or inst.empty:
        print(f"{EXPOSURE}: no instruments")
        return 2
    frames = []
    for key in PROTECTIVE + CONTRAST:
        r = ratios(inst, key)
        if r is None:
            print(f"  {key}: no cached outcome rows (run 01_ladder.py first)")
            continue
        frames.append(r)
    if len(frames) < 3:
        print("not enough endpoints")
        return 1
    w = pd.concat(frames, axis=1)
    print(f"{len(w)} shared variants across {w.shape[1]} endpoints")

    keys = list(w.columns)
    rows = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            s = w[[a, b]].dropna()
            if len(s) < 20:
                continue
            rows.append(dict(a=a, b=b, n=len(s),
                             r=float(np.corrcoef(s[a], s[b])[0, 1])))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "rung_concordance.csv", index=False, float_format="%.6g")
    print(f"-> {(OUT / 'rung_concordance.csv').relative_to(ROOT).as_posix()}\n")

    mat = pd.DataFrame(np.nan, index=keys, columns=keys)
    for _, r in df.iterrows():
        mat.loc[r["a"], r["b"]] = mat.loc[r["b"], r["a"]] = r["r"]
    for k in keys:
        mat.loc[k, k] = 1.0
    print(mat.round(2).to_string())

    prot = [k for k in PROTECTIVE if k in mat.index]
    con = [k for k in CONTRAST if k in mat.index]
    if len(prot) >= 2 and con:
        within = np.nanmean([mat.loc[a, b] for i, a in enumerate(prot)
                             for b in prot[i + 1:]])
        across = np.nanmean([mat.loc[a, b] for a in prot for b in con])
        print(f"\nmean r within the protective rungs: {within:.3f}")
        print(f"mean r between protective and null rungs: {across:.3f}")
        print("reading: shared per-variant effects across the protective rungs and lower "
              "agreement with MI / hypertension means one non-coronary axis, observed "
              "several times 鈥?not several independent findings")
    return 0


if __name__ == "__main__":
    sys.exit(main())

