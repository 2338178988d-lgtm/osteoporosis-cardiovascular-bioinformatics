"""Module 13A-diagnostics 鈥?why does DXA femoral neck BMD point the other way?

Script 01 returns a sign reversal: DXA femoral neck BMD -> FinnGen heart failure gives
b = +0.068 (p = 0.034) while heel-ultrasound eBMD gives b = -0.064 (p = 0.00057). A sign
reversal between two bone-density modalities is either a real dissociation or an artefact
of what the second instrument set happens to tag. This script separates those two.

Four checks:

  1. POSITIONAL CONCORDANCE. The two instrument sets share only 7 rsIDs, which says
     nothing, because fine-mapping picks different lead SNPs at the same locus. Match by
     position (+/- 500 kb) and correlate the two exposure effects at matched loci. A high
     correlation means the two files describe the same underlying signal and the reversal is
     not an allele-coding mistake.
  2. GENE CONTENT. Assign the DXA instruments to nearest genes (FinnGen `nearest_genes`
     column) and compare the gene list with the eBMD one. An adiposity/body-size-heavy list
     would explain a positive estimate.
  3. CONFOUNDING. MVMR of DXA BMD with BMI, lean mass and height. If the positive estimate
     collapses once body size is conditioned on, the DXA set is a body-size set.
  4. LOCUS REMOVAL. Drop DXA instruments within 500 kb of any BMI instrument and re-estimate.

Outputs (results/tables/13-robustness/):
  dxa_positional_concordance.csv
  dxa_instrument_genes.csv
  dxa_diagnostics_mvmr.csv

Usage:
  python analysis/13-robustness/scripts/04_dxa_diagnostics.py
  python analysis/13-robustness/scripts/04_dxa_diagnostics.py --self-test
"""
import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from _common import mr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / "results" / "tables" / "13-robustness"
FINNGEN_HF = ROOT / "data" / "raw" / "finngen" / "finngen_R12_I9_HEARTFAIL.gz"
WINDOW = 500_000


def align(b, ea, oa, ref_ea, ref_oa):
    if (ea, oa) == (ref_ea, ref_oa):
        return b
    if (ea, oa) == (ref_oa, ref_ea):
        return -b
    return None


def load(key):
    inst = mr.load_instruments(key)
    if inst is None:
        return None
    for c in ("chr", "pos", "bx", "sx", "eaf_x"):
        inst[c] = pd.to_numeric(inst[c], errors="coerce")
    return inst.dropna(subset=["chr", "pos"]).reset_index(drop=True)


def positional_join(a, b):
    """For each instrument in `a`, the closest instrument in `b` on the same chromosome."""
    rows = []
    for _, r in a.iterrows():
        same = b[b["chr"] == r["chr"]]
        if same.empty:
            continue
        d = (same["pos"] - r["pos"]).abs()
        j = d.idxmin()
        if d.loc[j] <= WINDOW:
            rows.append((r["SNP"], r["bx"], r["effect_allele"], r["other_allele"],
                         same.loc[j, "SNP"], same.loc[j, "bx"],
                         same.loc[j, "effect_allele"], same.loc[j, "other_allele"]))
    return pd.DataFrame(rows, columns=["snp_new", "bx_new", "ea_new", "oa_new",
                                       "snp_ref", "bx_ref", "ea_ref", "oa_ref"])


def scan_genes(snps, path=FINNGEN_HF):
    """Nearest gene for each rsID, read off the FinnGen heart-failure file."""
    want = set(snps)
    found = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        i_rs, i_gene = header.index("rsids"), header.index("nearest_genes")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= i_gene:
                continue
            if p[i_rs] in want and p[i_rs] not in found:
                found[p[i_rs]] = p[i_gene]
                if len(found) == len(want):
                    break
    return found


def mvmr(tab, extra):
    """MVMR-IVW of DXA BMD plus one extra exposure on the outcome."""
    sys.path.insert(0, str(ROOT / "analysis" / "05-hf-followup" / "scripts"))
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "fu", ROOT / "analysis" / "05-hf-followup" / "scripts" / "01_hf_followup.py")
    fu = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fu)
    cols = ["bx"] + [f"b_{k}" for k in extra]
    sub = tab.dropna(subset=cols + ["by", "sy"])
    if len(sub) < 20:
        return None
    X = sub[cols].to_numpy(float)
    res = fu.mvmr_ivw(X, sub["by"].to_numpy(float), sub["sy"].to_numpy(float))
    if res is None:
        return None
    beta, se, _ = res
    p = 2 * stats.norm.sf(np.abs(beta / se))
    return dict(nsnp=len(sub), beta=beta, se=se, p=p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--outcome", default="FG_HF")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    new, ref = load("FNBMD"), load("eBMD")
    if new is None or ref is None:
        print("need both instrument tables")
        return 2
    print(f"DXA instruments {len(new)}, eBMD instruments {len(ref)}")

    j = positional_join(new, ref)
    aligned = []
    for _, r in j.iterrows():
        b = align(r["bx_ref"], r["ea_ref"], r["oa_ref"], r["ea_new"], r["oa_new"])
        if b is not None:
            aligned.append((r["snp_new"], r["snp_ref"], r["bx_new"], b))
    ad = pd.DataFrame(aligned, columns=["snp_new", "snp_ref", "bx_new", "bx_ref_aligned"])
    print(f"positional matches within {WINDOW/1000:.0f} kb: {len(j)} of {len(new)} "
          f"({100*len(j)/len(new):.0f}%); {len(ad)} allele-alignable")
    if len(ad) >= 5:
        r = float(np.corrcoef(ad["bx_ref_aligned"], ad["bx_new"])[0, 1])
        slope = float(np.polyfit(ad["bx_ref_aligned"], ad["bx_new"], 1)[0])
        sign = float(np.mean(np.sign(ad["bx_ref_aligned"]) == np.sign(ad["bx_new"])))
        print(f"effect correlation at matched loci r={r:.3f}, slope={slope:.3f}, "
              f"sign agreement {sign:.0%}")
        pd.DataFrame([dict(window_kb=WINDOW/1000, n_dxa=len(new), n_matched=len(j),
                           n_alignable=len(ad), corr_matched=r, slope_matched=slope,
                           sign_agreement=sign)]).to_csv(
            OUT / "dxa_positional_concordance.csv", index=False, float_format="%.6g")
        ad.to_csv(OUT / "dxa_matched_loci.csv", index=False, float_format="%.6g")

    genes = scan_genes(new["SNP"].tolist())
    gdf = pd.DataFrame([(s, genes.get(s, "")) for s in new["SNP"]],
                       columns=["SNP", "nearest_gene"])
    gdf.to_csv(OUT / "dxa_instrument_genes.csv", index=False)
    print(f"genes assigned: {int((gdf.nearest_gene != '').sum())}/{len(gdf)}")
    print("  top genes:", ", ".join(gdf["nearest_gene"].head(18).tolist()))

    tab = None
    sys.path.insert(0, str(ROOT / "analysis" / "05-hf-followup" / "scripts"))
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "fu", ROOT / "analysis" / "05-hf-followup" / "scripts" / "01_hf_followup.py")
    fu = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fu)
    tab = fu.stack("FNBMD", a.outcome, ["BMI", "ALM", "HEIGHT"])
    if tab is None:
        print(f"{a.outcome}: no data")
        return 2
    base = mr.ivw(tab["bx"].to_numpy(float), tab["by"].to_numpy(float),
                  tab["sx"].to_numpy(float), tab["sy"].to_numpy(float), random=True)
    print(f"\nbaseline DXA BMD -> {a.outcome}: b={base['b']:+.4f} (p={base['p']:.3g}, "
          f"n={len(tab)})")
    rows = [dict(outcome=a.outcome, model="univariable", nsnp=len(tab), b=base["b"],
                 se=base["se"], p=base["p"], direction="+" if base["b"] > 0 else "-")]
    for extra in (["BMI"], ["ALM"], ["HEIGHT"], ["BMI", "HEIGHT"]):
        res = mvmr(tab, extra)
        if res is None:
            continue
        rows.append(dict(outcome=a.outcome, model="+".join(["DXA"] + extra), nsnp=res["nsnp"],
                         b=res["beta"][0], se=res["se"][0], p=res["p"][0],
                         direction="+" if res["beta"][0] > 0 else "-"))
        others = ", ".join(f"{k} {v:+.4f} (p={pv:.3g})"
                           for k, v, pv in zip(extra, res["beta"][1:], res["p"][1:]))
        print(f"  {'+'.join(extra):14s} n={res['nsnp']:3d}  DXA direct "
              f"{res['beta'][0]:+.4f} (p={res['p'][0]:.3g})   [{others}]")

    bmi_inst = mr.load_instruments("BMI")
    if bmi_inst is not None:
        bmi = bmi_inst.copy()
        bmi["chr"] = pd.to_numeric(bmi["chr"], errors="coerce")
        bmi["pos"] = pd.to_numeric(bmi["pos"], errors="coerce")
        pos = new.set_index("SNP")[["chr", "pos"]]
        tabp = tab.join(pos, on="SNP")
        keep = []
        for _, r in tabp.iterrows():
            near = bmi[(bmi["chr"] == r["chr"])
                       & ((bmi["pos"] - r["pos"]).abs() <= WINDOW)]
            keep.append(near.empty)
        keep = np.array(keep, bool)
        if keep.sum() >= 20:
            r2 = mr.ivw(tab["bx"].to_numpy(float)[keep], tab["by"].to_numpy(float)[keep],
                        tab["sx"].to_numpy(float)[keep], tab["sy"].to_numpy(float)[keep],
                        random=True)
            rows.append(dict(outcome=a.outcome, model="BMI-near loci removed",
                             nsnp=int(keep.sum()), b=r2["b"], se=r2["se"], p=r2["p"],
                             direction="+" if r2["b"] > 0 else "-"))
            print(f"  BMI-near loci removed ({int((~keep).sum())} dropped): n={int(keep.sum())} "
                  f"b={r2['b']:+.4f} (p={r2['p']:.3g})")
    pd.DataFrame(rows).to_csv(OUT / "dxa_diagnostics_mvmr.csv", index=False,
                              float_format="%.6g")
    print(f"\n-> {(OUT / 'dxa_diagnostics_mvmr.csv').relative_to(ROOT).as_posix()}")
    return 0


def selftest():
    """Allele alignment must flip the reference effect when the alleles are swapped."""
    v = align(0.05, "A", "G", "G", "A")
    w = align(0.05, "A", "G", "A", "G")
    u = align(0.05, "A", "G", "C", "T")
    print(f"selftest: swapped={v}, same={w}, incompatible={u}")
    ok = v == -0.05 and w == 0.05 and u is None
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

