"""FinnGen R12: independent-population validation of the eBMD -> heart-failure signal,
and a subtype test that can arbitrate the disagreement between the MR and NHANES layers.

Why FinnGen: eBMD and every GWAS used so far are dominated by UK Biobank, so no
existing "replication" was independent. FinnGen R12 is a Finnish cohort with no UK
Biobank overlap, which makes it the only genuinely independent test available without
a data application.

Three questions:
  1. REPLICATION. Do the 600 LD-clumped eBMD instruments reproduce the effect on
     `FG_HF` (heart failure, strict; 37,653 cases)?
  2. SUBTYPE. Does the effect survive in `FG_HF_CHD` (heart failure with coronary
     disease, 17,532 cases) and `FG_HF_OBESE` (heart failure with BMI >= 25, 22,965)?
     A BMD effect that is present in strict HF but absent from the coronary subtype is
     consistent with the MVMR finding that the effect is not routed through CAD; an
     effect confined to the obese subtype would contradict the BMI MVMR.
  3. REVERSE, ONCE. FinnGen HF as the exposure against the UK Biobank eBMD outcome.
     Both sides are then external to each other, which is the highest-powered clean
     version of the reverse direction. Run once and stop, whatever it shows.

Outputs (results/tables/05-hf-followup/):
  fg_replication.csv
  fg_subtype.csv
  fg_reverse.csv

Usage:
  python analysis/05-hf-followup/scripts/02_finngen.py
  python analysis/05-hf-followup/scripts/02_finngen.py --self-test
"""
import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
MR_DIR = HERE.parents[1] / "02-mr" / "scripts"
sys.path.insert(0, str(MR_DIR))

_spec = importlib.util.spec_from_file_location("mr", MR_DIR / "03_mr.py")
mr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mr)
import _sumstats

ROOT = HERE.parents[2]
OUT = ROOT / "results" / "tables" / "05-hf-followup"
EXPOSURE = "eBMD"


def run_pair(exp_key, out_key, min_snps=10):
    """Full set of estimators for one pair, using the module-02 machinery."""
    inst = mr.load_instruments(exp_key)
    if inst is None or inst.empty:
        return None, None
    out = mr.outcome_at_snps(out_key, inst["SNP"].tolist())
    if out is None or out.empty:
        return None, None
    d = mr.harmonise(inst, out)
    if len(d) < min_snps:
        print(f"  {exp_key} -> {out_key}: only {len(d)} usable SNPs")
        return None, None
    bx, by = d["bx"].to_numpy(float), d["by"].to_numpy(float)
    sx, sy = d["sx"].to_numpy(float), d["sy"].to_numpy(float)
    r = mr.ivw(bx, by, sx, sy, random=True)
    e = mr.egger(bx, by, sx, sy)
    row = dict(exposure=exp_key, outcome=out_key, nsnp=len(d),
               b=r["b"], se=r["se"], p=r["p"],
               lo=r["b"] - 1.96 * r["se"], hi=r["b"] + 1.96 * r["se"],
               or_=float(np.exp(r["b"])),
               or_lo=float(np.exp(r["b"] - 1.96 * r["se"])),
               or_hi=float(np.exp(r["b"] + 1.96 * r["se"])),
               q_p=r["q_p"],
               egger_b=e["b"], egger_se=e["se"], egger_p=e["p"],
               egger_intercept=e["intercept"],
               egger_intercept_p=e["intercept_p"])
    return row, d


def presso(d):
    """MR-PRESSO-style outlier removal, reusing module-02's implementation."""
    spec = importlib.util.spec_from_file_location("pl", MR_DIR / "04_pleiotropy.py")
    pl = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pl)
    bx, by = d["bx"].to_numpy(float), d["by"].to_numpy(float)
    sx, sy = d["sx"].to_numpy(float), d["sy"].to_numpy(float)
    keep, removed = pl.presso_outliers(bx, by, sx, sy)
    if keep.sum() < 10:
        return None, removed
    t = mr.ivw(bx[keep], by[keep], sx[keep], sy[keep], random=True)
    return t, removed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        print("no synthetic path here; the estimators are tested in 03_mr.py --self-test")
        return 0

    OUT.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------- 1. replication
    print("== 1. replication of eBMD -> HF in FinnGen (no UK Biobank overlap)")
    rows = []
    for key, label in [("HF", "UKB/HERMES (discovery, overlapping)"),
                       ("HF_SHAH", "Shah 2020 (overlapping)"),
                       ("FG_HF", "FinnGen R12 strict (independent)")]:
        r, _ = run_pair(EXPOSURE, key)
        if r is None:
            continue
        r["cohort"] = label
        rows.append(r)
        print(f"  {label:38s} n={r['nsnp']:4d}  b={r['b']:+.4f} "
              f"({r['lo']:+.4f}, {r['hi']:+.4f})  p={r['p']:.3g}  OR={r['or_']:.3f}  "
              f"Egger int p={r['egger_intercept_p']:.3g}")
    rep = pd.DataFrame(rows)
    rep.to_csv(OUT / "fg_replication.csv", index=False, float_format="%.6g")

    # trimmed estimate for the FinnGen pair as well, for a like-for-like comparison
    r_fg, d_fg = run_pair(EXPOSURE, "FG_HF")
    if d_fg is not None:
        t, removed = presso(d_fg)
        if t is not None:
            print(f"  FinnGen after removing {len(removed)} outlier variants: "
                  f"b={t['b']:+.4f} (p={t['p']:.3g}), n={len(d_fg) - len(removed)}")
            pd.DataFrame([dict(outcome="FG_HF", nsnp=len(d_fg), n_outliers=len(removed),
                               full_b=r_fg["b"], full_p=r_fg["p"],
                               trim_b=t["b"], trim_se=t["se"], trim_p=t["p"],
                               trim_or=float(np.exp(t["b"])))]) \
                .to_csv(OUT / "fg_trimmed.csv", index=False, float_format="%.6g")

    # ---------------------------------------------------------------- 2. subtype
    print("\n== 2. subtype test: is the effect present in the coronary or obese subtype?")
    rows = []
    for key, label in [("FG_HF", "all heart failure (37,653)"),
                       ("FG_HF_CHD", "HF + coronary disease (17,532)"),
                       ("FG_HF_OBESE", "HF + BMI>=25 (22,965)")]:
        r, _ = run_pair(EXPOSURE, key)
        if r is None:
            continue
        r["subtype"] = label
        rows.append(r)
        print(f"  {label:32s} n={r['nsnp']:4d}  b={r['b']:+.4f}  p={r['p']:.3g}  "
              f"OR={r['or_']:.3f}  q_p={r['q_p']:.3g}")
    sub = pd.DataFrame(rows)
    sub.to_csv(OUT / "fg_subtype.csv", index=False, float_format="%.6g")

    # ---------------------------------------------------------------- 3. reverse, once
    print("\n== 3. reverse direction, once: FinnGen HF -> UKB eBMD")
    reverse = []
    for key, label in [("FG_HF", "FinnGen HF -> eBMD")]:
        # build instruments from the FinnGen GWAS, then look them up in eBMD
        inst = build_finngen_instruments(key)
        if inst is None or inst.empty:
            print(f"  {label}: no instruments")
            continue
        out = mr.outcome_at_snps(EXPOSURE, inst["SNP"].tolist())
        if out is None or out.empty:
            print(f"  {label}: no outcome matches")
            continue
        d = mr.harmonise(inst, out)
        if len(d) < 10:
            print(f"  {label}: only {len(d)} usable SNPs")
            continue
        bx, by = d["bx"].to_numpy(float), d["by"].to_numpy(float)
        sx, sy = d["sx"].to_numpy(float), d["sy"].to_numpy(float)
        r = mr.ivw(bx, by, sx, sy, random=True)
        e = mr.egger(bx, by, sx, sy)
        reverse.append(dict(pair=label, nsnp=len(d), b=r["b"], se=r["se"], p=r["p"],
                            lo=r["b"] - 1.96 * r["se"], hi=r["b"] + 1.96 * r["se"],
                            q_p=r["q_p"], egger_intercept_p=e["intercept_p"]))
        print(f"  {label:26s} n={len(d):4d}  b={r['b']:+.4f} "
              f"({r['b'] - 1.96 * r['se']:+.4f}, {r['b'] + 1.96 * r['se']:+.4f})  p={r['p']:.3g}  "
              f"Egger int p={e['intercept_p']:.3g}")
    if reverse:
        pd.DataFrame(reverse).to_csv(OUT / "fg_reverse.csv", index=False, float_format="%.6g")
    return 0


def build_finngen_instruments(key, p_thresh=5e-8, kb=10000):
    """Extract genome-wide significant FinnGen variants, then clump with PLINK.

    Written here rather than reusing 02_extract_instruments.py because that script
    expects the GWAS Catalog column names and writes to the 02-mr table directory.
    """
    import subprocess
    import tempfile

    tab = OUT / f"inst_{key}_candidates.csv"
    if tab.exists():
        cand = pd.read_csv(tab)
    else:
        d = _sumstats.read_rows(key, p_thresh)
        if d.empty:
            return None
        for c in ["pos", "eaf", "beta", "se", "p", "n"]:
            d[c] = pd.to_numeric(d[c], errors="coerce")
        d = d.dropna(subset=["SNP", "beta", "se", "p", "eaf"])
        d = d[(d["eaf"] > 0.01) & (d["eaf"] < 0.99)]
        d = d[d["se"] > 0]
        d = d.sort_values("p").drop_duplicates("SNP")
        d = d[d["effect_allele"].str.len() == 1]
        d = d[d["other_allele"].str.len() == 1]
        d["F"] = (d["beta"] / d["se"]) ** 2
        d["r2"] = d["F"] / (d["F"] + d["n"] - 2)
        cand = d.reset_index(drop=True)
        cand.to_csv(tab, index=False)
    print(f"  {key}: {len(cand)} candidates at p<{p_thresh:g}")

    ref = ROOT / "data" / "raw" / "reference"
    plink = ref / "plink" / "plink.exe"
    if not plink.exists() or not (ref / "EUR.bed").exists():
        print("  PLINK or EUR reference missing; using distance pruning")
        return prune_distance(cand, kb)

    with tempfile.TemporaryDirectory() as tmp:
        assoc = Path(tmp) / f"{key}.assoc"
        with open(assoc, "w", encoding="utf-8") as f:
            f.write("CHR SNP BP P\n")
            for _, r in cand.iterrows():
                f.write(f"{int(r['chr'])} {r['SNP']} {int(r['pos'])} {r['p']:.6g}\n")
        out = Path(tmp) / f"{key}_c"
        cmd = [str(plink), "--bfile", str(ref / "EUR"), "--clump", str(assoc),
               "--clump-p1", "1", "--clump-r2", "0.001", "--clump-kb", str(kb),
               "--clump-snp-field", "SNP", "--clump-field", "P",
               "--out", str(out), "--allow-no-sex", "--memory", "4000"]
        subprocess.run(cmd, capture_output=True, text=True)
        clumped = Path(f"{out}.clumped")
        if not clumped.exists():
            return prune_distance(cand, kb)
        c = pd.read_csv(clumped, sep=r"\s+", engine="python")
        keep = set(c["SNP"].astype(str))
        inst = cand[cand["SNP"].isin(keep)].sort_values("p").reset_index(drop=True)
    inst = inst.rename(columns={"beta": "bx", "se": "sx", "p": "px", "eaf": "eaf_x", "n": "n_x"})
    inst.to_csv(OUT / f"instruments_{key}.csv", index=False)
    print(f"  {key}: {len(inst)} LD-independent instruments, min F = {inst['F'].min():.1f}")
    return inst


def prune_distance(cand, kb):
    kept, dropped = [], {}
    for chrom, g in cand.sort_values("p").groupby("chr", sort=False):
        leads = []
        for _, row in g.iterrows():
            hit = next((l for l in leads if abs(l[0] - row["pos"]) <= kb * 1000), None)
            if hit is None:
                leads.append((row["pos"], row["SNP"]))
                kept.append(row)
            else:
                dropped[row["SNP"]] = hit[1]
    d = pd.DataFrame(kept).sort_values("p").reset_index(drop=True)
    return d.rename(columns={"beta": "bx", "se": "sx", "p": "px", "eaf": "eaf_x", "n": "n_x"})


if __name__ == "__main__":
    sys.exit(main())

