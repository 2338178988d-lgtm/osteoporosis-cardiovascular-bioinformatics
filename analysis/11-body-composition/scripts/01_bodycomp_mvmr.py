"""Module 11 鈥?is eBMD a bone signal or a body-composition signal?

eBMD is a DXA image phenotype: it is bone mineral density measured in a limb, and its
genetic architecture is correlated with lean mass and with stature. That leaves one
alternative explanation for eBMD -> heart failure that nothing in this project has tested:

    variants that raise eBMD partly raise muscle mass and body size; low muscle mass
    (sarcopenia) raises heart-failure risk through frailty, reduced cardiac reserve and
    sarcopenic obesity. The apparent bone -> heart axis would then be a muscle -> heart
    axis seen through a bone-shaped instrument.

This module separates the two by multivariable MR against appendicular lean mass (ALM,
n=450,243) and standing height (n=673,878), both UK Biobank GWAS on the same harmonised
platform as the eBMD exposure, so the instruments and LD structure are comparable.

Four models per outcome:
  M1  eBMD alone                       (baseline, reproduces modules 02/05)
  M2  eBMD + ALM                       (is the eBMD effect muscle?)
  M3  eBMD + height                    (is it body size?)
  M4  eBMD + ALM + BMI                 (does lean mass survive adjustment for total mass?)

plus the univariable ALM -> HF estimate on ALM's own instruments: does muscle itself move
heart failure, independently of bone?

Decision rule, fixed before the run:
  eBMD direct effect attenuates  < 10 %   the bone reading stands
                                 10-30 %  bone and muscle share the axis
                                 > 30 %   eBMD is a proxy and muscle is the exposure

Outputs (results/tables/11-body-composition/):
  bodycomp_mvmr.csv        direct effect per model, per outcome
  bodycomp_diagnostics.csv conditional F and exposure-exposure correlation per model
  bodycomp_alm_own.csv     univariable ALM -> HF on ALM's own instruments

Usage:
  python analysis/11-body-composition/scripts/01_bodycomp_mvmr.py
  python analysis/11-body-composition/scripts/01_bodycomp_mvmr.py --self-test
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
for name, path in [("mr", MR_DIR / "03_mr.py"),
                   ("fu", HERE.parents[1] / "05-hf-followup" / "scripts" / "01_hf_followup.py")]:
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    globals()[name] = mod

ROOT = HERE.parents[2]
OUT = ROOT / "results" / "tables" / "11-body-composition"
EXPOSURE = "eBMD"
OUTCOMES = {"HF": "HERMES (overlapping)", "FG_HF": "FinnGen R12 (independent)"}
MODELS = [("M1_eBMD", []),
          ("M2_eBMD_ALM", ["ALM"]),
          ("M3_eBMD_HEIGHT", ["HEIGHT"]),
          ("M4_eBMD_ALM_BMI", ["ALM", "BMI"])]


def fit(tab, extra):
    """MVMR-IVW of eBMD plus the extra exposures on the outcome; eBMD coefficient first."""
    cols = ["bx"] + [f"b_{k}" for k in extra]
    sub = tab.dropna(subset=cols + ["by", "sy"])
    if len(sub) < 20:
        return None, sub
    X = sub[cols].to_numpy(float)
    y = sub["by"].to_numpy(float)
    sy = sub["sy"].to_numpy(float)
    res = fu.mvmr_ivw(X, y, sy)
    if res is None:
        return None, sub
    beta, se, _ = res
    p = 2 * stats.norm.sf(np.abs(beta / se))
    return dict(nsnp=len(sub), beta=beta, se=se, p=p), sub


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return selftest()

    OUT.mkdir(parents=True, exist_ok=True)
    rows, diag, own = [], [], []
    for key, label in OUTCOMES.items():
        print(f"\n=== {key} ({label})")
        tab = fu.stack(EXPOSURE, key, ["ALM", "HEIGHT", "BMI"])
        if tab is None:
            print("  no data")
            continue
        base = mr.ivw(tab["bx"].to_numpy(float), tab["by"].to_numpy(float),
                      tab["sx"].to_numpy(float), tab["sy"].to_numpy(float), random=True)
        print(f"  baseline eBMD -> {key}: b={base['b']:+.4f} (p={base['p']:.3g}, "
              f"n={len(tab)})")
        for name, extra in MODELS:
            if not extra:
                rows.append(dict(outcome=key, model=name, exposures="eBMD", nsnp=len(tab),
                                 direct_b=base["b"], direct_se=base["se"], direct_p=base["p"],
                                 attenuation_pct=0.0))
                continue
            res, sub = fit(tab, extra)
            if res is None:
                print(f"  {name}: too few shared instruments")
                continue
            att = 100 * (1 - res["beta"][0] / base["b"]) if base["b"] else np.nan
            rows.append(dict(outcome=key, model=name, exposures="+".join(["eBMD"] + extra),
                             nsnp=res["nsnp"], direct_b=res["beta"][0], direct_se=res["se"][0],
                             direct_p=res["p"][0], attenuation_pct=att))
            # collinearity: conditional F for eBMD given the extras, and the raw correlation
            cf = fu.conditional_f(sub["bx"].to_numpy(float), sub["sx"].to_numpy(float),
                                  sub[f"b_{extra[0]}"].to_numpy(float), len(sub))
            corr = np.corrcoef(sub["bx"].to_numpy(float),
                               sub[f"b_{extra[0]}"].to_numpy(float))[0, 1]
            diag.append(dict(outcome=key, model=name, nsnp=res["nsnp"], conditional_F_eBMD=cf,
                             corr_eBMD_first_extra=float(corr),
                             extra_betas="; ".join(f"{k}:{v:+.4f}" for k, v in
                                                   zip(extra, res["beta"][1:])),
                             extra_p="; ".join(f"{k}:{v:.3g}" for k, v in
                                               zip(extra, res["p"][1:]))))
            print(f"  {name:16s} n={res['nsnp']:4d}  direct {res['beta'][0]:+.4f} "
                  f"(p={res['p'][0]:.3g})  attenuation {att:+6.1f}%  "
                  f"conditional F {cf:.0f}")

        # does muscle move heart failure on its own instruments, and is that just body size?
        tab_alm = fu.stack("ALM", key, ["BMI", "HEIGHT"])
        if tab_alm is not None and len(tab_alm) >= 20:
            uni = mr.ivw(tab_alm["bx"].to_numpy(float), tab_alm["by"].to_numpy(float),
                         tab_alm["sx"].to_numpy(float), tab_alm["sy"].to_numpy(float),
                         random=True)
            own.append(dict(outcome=key, exposure="ALM", model="univariable",
                            nsnp=len(tab_alm), b=uni["b"], se=uni["se"], p=uni["p"],
                            or_=float(np.exp(uni["b"])), attenuation_pct=0.0))
            print(f"  ALM own instruments -> {key}: n={len(tab_alm)} b={uni['b']:+.4f} "
                  f"(p={uni['p']:.3g})")
            for extra in (["BMI"], ["HEIGHT"], ["BMI", "HEIGHT"]):
                res, sub = fit(tab_alm, extra)
                if res is None:
                    continue
                att = 100 * (1 - res["beta"][0] / uni["b"]) if uni["b"] else np.nan
                own.append(dict(outcome=key, exposure="ALM", model="+".join(["ALM"] + extra),
                                nsnp=res["nsnp"], b=res["beta"][0], se=res["se"][0],
                                p=res["p"][0], or_=float(np.exp(res["beta"][0])),
                                attenuation_pct=att))
                print(f"    adjusted for {'+'.join(extra):12s} n={res['nsnp']:4d} "
                      f"ALM direct {res['beta'][0]:+.4f} (p={res['p'][0]:.3g})  "
                      f"attenuation {att:+6.1f}%")
    if rows:
        pd.DataFrame(rows).to_csv(OUT / "bodycomp_mvmr.csv", index=False, float_format="%.6g")
        print(f"\n-> {(OUT / 'bodycomp_mvmr.csv').relative_to(ROOT).as_posix()}")
    if diag:
        pd.DataFrame(diag).to_csv(OUT / "bodycomp_diagnostics.csv", index=False,
                                  float_format="%.6g")
    if own:
        pd.DataFrame(own).to_csv(OUT / "bodycomp_alm_own.csv", index=False,
                                 float_format="%.6g")
    return 0


def selftest():
    """A mediator that truly carries part of the effect must show intermediate attenuation.

    Simulated instruments move eBMD and ALM together (conditional F 鈮?15, comparable to
    the real models); the outcome is driven by both. The estimator has to recover the
    planted direct effect and land in the "shared axis" band of the decision rule.
    """
    rng = np.random.default_rng(7)
    n = 4000
    bx = rng.normal(0.05, 0.02, n)
    bm = 0.3 * bx + rng.normal(0, 0.03, n)
    sy = np.full(n, 0.01)
    by = -0.02 * bx - 0.03 * bm + rng.normal(0, sy)
    beta, se, _ = fu.mvmr_ivw(np.column_stack([bx, bm]), by, sy)
    p = 2 * stats.norm.sf(np.abs(beta / se))
    cf = fu.conditional_f(bx, np.full(n, 0.005), bm, n)
    uni = mr.ivw(bx, by, np.full(n, 0.005), sy, random=True)
    att = 100 * (1 - beta[0] / uni["b"])
    print(f"selftest: planted direct -0.020, estimated {beta[0]:+.4f} (p={p[0]:.1e}); "
          f"univariable {uni['b']:+.4f}; attenuation {att:+.1f}%; conditional F {cf:.0f}")
    ok = abs(beta[0] + 0.02) < 0.006 and p[0] < 0.01 and 15 < att < 60
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

