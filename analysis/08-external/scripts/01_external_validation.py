"""Module 08 鈥?external validation in a third population, and by ejection fraction.

Two questions module 05 could not answer.

1. THIRD POPULATION. The independent replication so far is FinnGen, a Finnish founder
   cohort, where the effect was *larger* than in UK Biobank (OR 0.938 vs 0.967). That
   raises the possibility of a population-specific effect. BioBank Japan provides a cohort
   with no UK Biobank and no Finnish overlap:
     FG_EA_HF  chronic heart failure, East Asian, GCST90018586 (Sakaue 2021)
     FG_EA_OP  osteoporosis, East Asian, GCST90018667

2. EJECTION FRACTION. Module 05's subtype test was by comorbidity (coronary disease,
   obesity), because FinnGen publishes no ejection-fraction split. HFpEF and HFrEF GWAS
   now exist in two populations:
     HFPEUR / HFREUR  non-ischaemic HFpEF / HFrEF, European (Henry 2025)
     HFPEA / HFREA    HFpEF / HFrEF, East Asian (Enzan 2025)
   A BMD effect present in HFpEF but not HFrEF would be a mechanistic prediction that the
   obesity-adjusted MVMR result does not make; an effect in both would say ejection
   fraction is not the axis.

Outputs (results/tables/08-external/):
  ea_replication.csv    eBMD -> HF in UKB, FinnGen and BioBank Japan
  ea_heterogeneity.csv  Cochran's Q across the three cohorts
  ea_reverse.csv        BioBank Japan HF -> UKB eBMD
  ef_subtype.csv        eBMD -> HFpEF and HFrEF in both populations
  ef_contrast.csv       HFpEF vs HFrEF contrast within each population

Usage:
  python analysis/08-external/scripts/01_external_validation.py
  python analysis/08-external/scripts/01_external_validation.py --self-test
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

ROOT = HERE.parents[2]
OUT = ROOT / "results" / "tables" / "08-external"
EXPOSURE = "eBMD"


def have(key):
    """True when a trait's summary-statistics file is present on disk."""
    import _sumstats
    return _sumstats.FILES.get(key) is not None and _sumstats.path_of(key).exists()


def estimate(exp_key, out_key, inst=None, label=""):
    """IVW random-effects estimate for one pair; returns a row dict or None."""
    if not have(out_key):
        print(f"  {label or out_key}: file not downloaded, skipped")
        return None
    if inst is None:
        inst = mr.load_instruments(exp_key)
    if inst is None or inst.empty:
        print(f"  {label or exp_key}: no instruments")
        return None
    out = mr.outcome_at_snps(out_key, inst["SNP"].tolist())
    if out is None or out.empty:
        print(f"  {label or out_key}: no outcome matches")
        return None
    d = mr.harmonise(inst, out)
    if len(d) < 10:
        print(f"  {label or out_key}: only {len(d)} usable SNPs")
        return None
    bx, by = d["bx"].to_numpy(float), d["by"].to_numpy(float)
    sx, sy = d["sx"].to_numpy(float), d["sy"].to_numpy(float)
    r = mr.ivw(bx, by, sx, sy, random=True)
    e = mr.egger(bx, by, sx, sy)
    row = dict(pair=label or f"{exp_key} -> {out_key}", exposure=exp_key, outcome=out_key,
               nsnp=len(d), b=r["b"], se=r["se"], p=r["p"],
               lo=r["b"] - 1.96 * r["se"], hi=r["b"] + 1.96 * r["se"],
               or_=float(np.exp(r["b"])),
               or_lo=float(np.exp(r["b"] - 1.96 * r["se"])),
               or_hi=float(np.exp(r["b"] + 1.96 * r["se"])),
               q_p=r["q_p"], egger_intercept_p=e["intercept_p"])
    print(f"  {row['pair']:38s} n={len(d):4d}  b={r['b']:+.4f} "
          f"({row['lo']:+.4f}, {row['hi']:+.4f})  p={r['p']:.3g}  OR={row['or_']:.3f}  "
          f"Egger int p={e['intercept_p']:.3g}")
    return row


def external_instruments(key):
    """Instrument set for a FinnGen / BioBank Japan trait, via the module-05 helper."""
    spec = importlib.util.spec_from_file_location(
        "fgen", HERE.parents[1] / "05-hf-followup" / "scripts" / "02_finngen.py")
    fgen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fgen)
    return fgen.build_finngen_instruments(key)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    inst = mr.load_instruments(EXPOSURE)

    print("== 1. third population: eBMD -> heart failure across three cohorts")
    rows = []
    for key, label in [("HF", "UKB / HERMES (overlapping with exposure)"),
                       ("FG_HF", "FinnGen R12 (Finnish)"),
                       ("FG_EA_HF", "BioBank Japan (East Asian)")]:
        r = estimate(EXPOSURE, key, inst=inst, label=label)
        if r:
            rows.append(r)
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(OUT / "ea_replication.csv", index=False, float_format="%.6g")
        if len(df) >= 2:
            w = 1 / df["se"] ** 2
            pooled = float((w * df["b"]).sum() / w.sum())
            q = float((w * (df["b"] - pooled) ** 2).sum())
            dfq = len(df) - 1
            p_q = float(stats.chi2.sf(q, dfq))
            print(f"  heterogeneity across cohorts: Q={q:.2f}, df={dfq}, p={p_q:.3g}")
            pd.DataFrame([dict(n_cohorts=len(df), pooled_b=pooled, Q=q, df=dfq,
                               p_heterogeneity=p_q)]).to_csv(
                OUT / "ea_heterogeneity.csv", index=False, float_format="%.6g")

    print("\n== 2. reverse, East Asian population")
    if have("FG_EA_HF"):
        inst_ea = external_instruments("FG_EA_HF")
        if inst_ea is not None and not inst_ea.empty:
            r = estimate("FG_EA_HF", EXPOSURE, inst=inst_ea, label="BBJ HF -> UKB eBMD")
            if r:
                pd.DataFrame([r]).to_csv(OUT / "ea_reverse.csv", index=False,
                                         float_format="%.6g")
    else:
        print("  FG_EA_HF not downloaded, skipped")

    print("\n== 3. ejection fraction: eBMD -> HFpEF vs HFrEF")
    ef = []
    for key, label in [("HFPEUR", "HFpEF European (3,590 cases)"),
                       ("HFREUR", "HFrEF European (4,975 cases)"),
                       ("HFPEA", "HFpEF East Asian"),
                       ("HFREA", "HFrEF East Asian")]:
        r = estimate(EXPOSURE, key, inst=inst, label=label)
        if r:
            ef.append(r)
    if ef:
        dfe = pd.DataFrame(ef)
        dfe.to_csv(OUT / "ef_subtype.csv", index=False, float_format="%.6g")
        contrasts = []
        for pop, p_key, r_key in [("European", "HFPEUR", "HFREUR"),
                                  ("East Asian", "HFPEA", "HFREA")]:
            s = dfe[dfe.outcome == p_key]
            t = dfe[dfe.outcome == r_key]
            if s.empty or t.empty:
                continue
            s, t = s.iloc[0], t.iloc[0]
            diff = s["b"] - t["b"]
            # conservative: treat the two outcome estimates as independent, which
            # over-states the SE because they share the same instruments
            se = float(np.hypot(s["se"], t["se"]))
            z = diff / se
            contrasts.append(dict(population=pop, b_hfpef=s["b"], b_hfref=t["b"],
                                  diff=diff, se_upper=se, z=z,
                                  p_upper=float(2 * stats.norm.sf(abs(z)))))
            print(f"  {pop}: HFpEF {s['b']:+.4f} vs HFrEF {t['b']:+.4f}, "
                  f"diff {diff:+.4f} (p<={2 * stats.norm.sf(abs(z)):.3g})")
        if contrasts:
            pd.DataFrame(contrasts).to_csv(OUT / "ef_contrast.csv", index=False,
                                           float_format="%.6g")
    return 0


def selftest():
    """The cohort-heterogeneity test must be ~0 for identical cohorts and large for
    divergent ones."""
    def qtest(bs, se=0.02):
        w = np.full(len(bs), 1 / se ** 2)
        b = np.asarray(bs)
        pooled = (w * b).sum() / w.sum()
        q = float((w * (b - pooled) ** 2).sum())
        return q, float(stats.chi2.sf(q, len(b) - 1))
    q_same, p_same = qtest([-0.04, -0.04, -0.04])
    q_diff, p_diff = qtest([-0.04, -0.30, 0.05])
    print(f"selftest: identical cohorts Q={q_same:.2f} (p={p_same:.2f})")
    print(f"selftest: divergent cohorts Q={q_diff:.1f} (p={p_diff:.2g})")
    ok = q_same < 0.01 and p_same > 0.9 and q_diff > 50 and p_diff < 1e-6
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

