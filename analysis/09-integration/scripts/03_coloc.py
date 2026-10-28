"""Module 09 step 3 鈥?colocalisation: is the bone-density signal and the heart-failure
signal at the same causal variant?

A shared nearest gene is not evidence of a shared mechanism. Two traits can associate at
the same locus because of the same variant (colocalisation), because of two distinct
variants in linkage disequilibrium, or for unrelated reasons in the same genomic region.
Only the first supports a mechanism. This module computes the posterior probability that
eBMD and heart failure share a causal variant at each of the instrument loci.

Method: coloc-abf (Giambartolomei 2014) with at most one causal variant per trait per
locus. The Wakefield approximate Bayes factor for each SNP is

    z = beta / se,  V = se^2,  W = (sd_prior * sd_Y)^2
    log ABF = 0.5 * ( log(V/(V+W)) + z^2 * W/(V+W) )

with sd_prior = 0.15 for a quantitative trait and 0.2 for a case-control trait, and
sd_Y = 1/sqrt(s(1-s)) for a case-control trait with case fraction s. Priors are the coloc
defaults p1 = p2 = 1e-4, p12 = 1e-5. The four hypotheses are combined in log space with
the (n-1) and (n-2) null factors retained, so PP.H4 is not an artefact of a shortcut.

Only the 600 instrument loci are examined, and only +/- 500 kb around each. Both files are
scanned once per chromosome set rather than once per locus.

Outputs (results/tables/09-integration/):
  coloc_loci.csv      per locus: PP.H0..H4, n SNPs, top SNP, whether H4 is supported
  coloc_summary.csv   counts by posterior class and the driver loci

Usage:
  python analysis/09-integration/scripts/03_coloc.py
  python analysis/09-integration/scripts/03_coloc.py --self-test
"""
import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parents[3]
EBMD = ROOT / "data" / "raw" / "gwas" / "GCST90029004_eBMD.h.tsv.gz"
FG_HF = ROOT / "data" / "raw" / "finngen" / "finngen_R12_I9_HEARTFAIL.gz"
INSTR = ROOT / "results" / "tables" / "02-mr" / "instruments_eBMD.csv"
OUT = ROOT / "results" / "tables" / "09-integration"

WINDOW = 500_000
N_CASES, N_CONTROLS = 37_653, 462_695
P1 = P2 = 1e-4
P12 = 1e-5
LOG10 = np.log(10)


def sd_y_casecontrol(cases, controls):
    s = cases / (cases + controls)
    return 1.0 / np.sqrt(s * (1 - s))


def log_abf(beta, se, w):
    """Wakefield log10 ABF; w is the prior variance of the true effect."""
    v = se ** 2
    z = beta / se
    return 0.5 * (np.log(v / (v + w)) + z ** 2 * w / (v + w)) / LOG10


def coloc_abf(labf1, labf2, n, p1=P1, p2=P2, p12=P12):
    """Posterior probabilities of H0..H4; single causal variant per trait."""
    # convert to natural log for stable summation, then back
    a1 = labf1 * LOG10
    a2 = labf2 * LOG10
    s1 = logsumexp(a1)
    s2 = logsumexp(a2)
    s12 = logsumexp(a1 + a2)
    # H3 term: (sum ABF1)(sum ABF2) - sum(ABF1*ABF2), computed in log space
    both = s1 + s2
    diff = np.log(np.expm1(both - s12)) + s12 if both > s12 else -np.inf
    lh = np.array([
        0.0,                                                     # H0
        np.log(p1) + s1 + (n - 1) * np.log1p(-p1) + n * np.log1p(-p2),
        np.log(p2) + s2 + n * np.log1p(-p1) + (n - 1) * np.log1p(-p2),
        np.log(p1) + np.log(p2) + diff + (n - 2) * np.log1p(-p1) + (n - 2) * np.log1p(-p2),
        np.log(p12) + s12 + (n - 1) * np.log1p(-p1) + (n - 1) * np.log1p(-p2),
    ])
    post = np.exp(lh - logsumexp(lh))
    return post


def build_windows(inst, window=WINDOW):
    """Merge instrument loci into non-overlapping windows per chromosome."""
    w = inst[["chr", "pos"]].copy()
    w["chr"] = w["chr"].astype(int)
    w["start"] = (w["pos"] - window).clip(lower=1)
    w["end"] = w["pos"] + window
    out = []
    for c, g in w.groupby("chr"):
        g = g.sort_values("start")
        cur_s, cur_e = int(g.iloc[0]["start"]), int(g.iloc[0]["end"])
        for _, r in g.iloc[1:].iterrows():
            if r["start"] <= cur_e:
                cur_e = max(cur_e, int(r["end"]))
            else:
                out.append((c, cur_s, cur_e))
                cur_s, cur_e = int(r["start"]), int(r["end"])
        out.append((c, cur_s, cur_e))
    return out


def in_windows(chrom, pos, lookup, wins):
    wins_chr = lookup.get(int(chrom))
    if wins_chr is None:
        return False
    for s, e in wins_chr:
        if s <= pos <= e:
            return True
    return False


def scan_ebmd(wins):
    lookup = {}
    for c, s, e in wins:
        lookup.setdefault(c, []).append((s, e))
    rows = []
    for chunk in pd.read_csv(EBMD, sep="\t", usecols=["hm_rsid", "hm_chrom", "hm_pos",
                                                      "hm_effect_allele", "hm_other_allele",
                                                      "hm_effect_allele_frequency",
                                                      "hm_beta", "standard_error", "p_value", "n"],
                             chunksize=2_000_000, low_memory=False):
        for c, g in chunk.groupby("hm_chrom"):
            try:
                g2 = g[g["hm_pos"].astype(float).map(lambda p: in_windows(c, p, lookup, wins))]
            except (ValueError, TypeError):
                continue
            if len(g2):
                rows.append(g2)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def scan_fg(wins):
    lookup = {}
    for c, s, e in wins:
        lookup.setdefault(c, []).append((s, e))
    rows = []
    for chunk in pd.read_csv(FG_HF, sep="\t", usecols=["#chrom", "pos", "ref", "alt", "rsids",
                                                       "pval", "beta", "sebeta", "af_alt"],
                             chunksize=2_000_000, low_memory=False):
        for c, g in chunk.groupby("#chrom"):
            try:
                g2 = g[g["pos"].astype(float).map(lambda p: in_windows(c, p, lookup, wins))]
            except (ValueError, TypeError):
                continue
            if len(g2):
                rows.append(g2)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--min-snps", type=int, default=50)
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    inst = pd.read_csv(INSTR)
    wins = build_windows(inst)
    print(f"{len(inst)} instruments -> {len(wins)} merged windows of +/-{WINDOW // 1000} kb")

    print("scanning eBMD ...")
    e = scan_ebmd(wins)
    print(f"  {len(e):,} variants in windows")
    print("scanning FinnGen HF ...")
    f = scan_fg(wins)
    print(f"  {len(f):,} variants in windows")

    e = e.rename(columns={"hm_rsid": "SNP", "hm_chrom": "chr", "hm_pos": "pos",
                          "hm_effect_allele": "ea", "hm_other_allele": "oa",
                          "hm_effect_allele_frequency": "eaf", "hm_beta": "b1",
                          "standard_error": "s1", "p_value": "p1", "n": "n1"})
    f = f.rename(columns={"#chrom": "chr", "pos": "pos", "ref": "oa2", "alt": "ea2",
                          "rsids": "SNP", "pval": "p2", "beta": "b2", "sebeta": "s2",
                          "af_alt": "eaf2"})
    for d in (e, f):
        d["chr"] = pd.to_numeric(d["chr"], errors="coerce")
        d["pos"] = pd.to_numeric(d["pos"], errors="coerce")
    e = e.dropna(subset=["chr", "pos", "b1", "s1"])
    f = f.dropna(subset=["chr", "pos", "b2", "s2"])

    m = e.merge(f.drop(columns=["chr", "pos"]), on="SNP", how="inner")
    print(f"  {len(m):,} variants shared by both traits")
    # align the outcome to the exposure effect allele
    same = (m["ea"] == m["ea2"]) & (m["oa"] == m["oa2"])
    flip = (m["ea"] == m["oa2"]) & (m["oa"] == m["ea2"])
    keep = same | flip
    m = m[keep].copy()
    m.loc[flip[keep], "b2"] *= -1
    m.loc[flip[keep], "eaf2"] = 1 - m.loc[flip[keep], "eaf2"]
    pal = m.apply(lambda r: {r["ea"], r["oa"]} in ({"A", "T"}, {"C", "G"}), axis=1)
    amb = pal & (m["eaf2"].between(0.42, 0.58))
    m = m[~amb]
    print(f"  {len(m):,} after allele alignment and palindromic filter")

    w1 = (0.15) ** 2
    w2 = (0.2 * sd_y_casecontrol(N_CASES, N_CONTROLS)) ** 2
    m["labf1"] = log_abf(m["b1"].to_numpy(float), m["s1"].to_numpy(float), w1)
    m["labf2"] = log_abf(m["b2"].to_numpy(float), m["s2"].to_numpy(float), w2)

    rows = []
    for _, r in inst.iterrows():
        c, p = int(r["chr"]), int(r["pos"])
        win = m[(m["chr"].astype(int) == c) & (m["pos"].between(p - WINDOW, p + WINDOW))]
        if len(win) < a.min_snps:
            rows.append(dict(SNP=r["SNP"], chr=c, pos=p, n_snps=len(win),
                             PP_H0=np.nan, PP_H1=np.nan, PP_H2=np.nan, PP_H3=np.nan,
                             PP_H4=np.nan, top_snp="", note="too few shared variants"))
            continue
        post = coloc_abf(win["labf1"].to_numpy(), win["labf2"].to_numpy(), len(win))
        top = win.loc[(win["labf1"] + win["labf2"]).idxmax()]
        rows.append(dict(SNP=r["SNP"], chr=c, pos=p, n_snps=len(win),
                         PP_H0=post[0], PP_H1=post[1], PP_H2=post[2],
                         PP_H3=post[3], PP_H4=post[4],
                         top_snp=top["SNP"], note=""))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "coloc_loci.csv", index=False, float_format="%.6g")

    ok = d.dropna(subset=["PP_H4"])
    print(f"\nloci with a usable colocalisation estimate: {len(ok)}/{len(d)}")
    for thr in (0.5, 0.8):
        n4 = int((ok["PP_H4"] >= thr).sum())
        n3 = int((ok["PP_H3"] >= thr).sum())
        print(f"  PP.H4 >= {thr}: {n4} loci | PP.H3 >= {thr}: {n3} loci")

    summ = pd.DataFrame([
        dict(metric="loci_tested", value=len(d)),
        dict(metric="loci_usable", value=len(ok)),
        dict(metric="PP_H4_ge_0.8", value=int((ok["PP_H4"] >= 0.8).sum())),
        dict(metric="PP_H4_ge_0.5", value=int((ok["PP_H4"] >= 0.5).sum())),
        dict(metric="PP_H3_ge_0.8", value=int((ok["PP_H3"] >= 0.8).sum())),
        dict(metric="median_PP_H0", value=float(ok["PP_H0"].median())),
        dict(metric="median_PP_H4", value=float(ok["PP_H4"].median())),
    ])
    summ.to_csv(OUT / "coloc_summary.csv", index=False, float_format="%.6g")
    print("\nsummary:")
    print(summ.to_string(index=False))

    print("\ntop colocalising loci (PP.H4):")
    print(ok.nlargest(12, "PP_H4")[["SNP", "chr", "pos", "n_snps", "PP_H3", "PP_H4", "top_snp"]]
          .round(4).to_string(index=False))
    return 0


def selftest():
    """A shared causal variant must give PP.H4 ~ 1; two distinct causal variants in the
    same window must give PP.H3 ~ 1."""
    rng = np.random.default_rng(11)
    n = 400
    pos = np.arange(n) * 1000
    # shared signal: the SAME variant drives both traits
    i = 200
    b1 = rng.normal(0, 0.01, n)
    b2 = rng.normal(0, 0.01, n)
    b1[i], b2[i] = 0.18, 0.16
    s1 = np.full(n, 0.02)
    s2 = np.full(n, 0.02)
    w1 = 0.15 ** 2
    w2 = (0.2 * sd_y_casecontrol(N_CASES, N_CONTROLS)) ** 2
    p_shared = coloc_abf(log_abf(b1, s1, w1), log_abf(b2, s2, w2), n)

    # distinct signals: different variants
    b1b = rng.normal(0, 0.01, n)
    b2b = rng.normal(0, 0.01, n)
    b1b[80] = 0.18
    b2b[320] = 0.16
    p_distinct = coloc_abf(log_abf(b1b, s1, w1), log_abf(b2b, s2, w2), n)

    print(f"selftest shared:    H3={p_shared[3]:.3f} H4={p_shared[4]:.3f}")
    print(f"selftest distinct:  H3={p_distinct[3]:.3f} H4={p_distinct[4]:.3f}")
    ok = p_shared[4] > 0.8 and p_distinct[3] > 0.8
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

