"""Select genetic instruments for the Mendelian-randomisation layer (02-mr).

Reads the GWAS Catalog harmonised files in data/raw/gwas/ (column names differ between
submitting studies, so they are normalised by _sumstats.py) and writes one instrument
table per trait:

    results/tables/02-mr/instruments_<KEY>.csv

Steps: keep genome-wide significant variants, drop low-MAF and malformed rows, then make
the instruments independent either by 1000G EUR LD clumping (Ensembl REST, r2 < 0.001
within 500 kb) or by positional pruning. LD clumping needs one HTTP request per SNP pair,
so it is preceded by a 100 kb distance pre-prune; a full LD pass over the raw candidate
list is not affordable this way and needs a local 1000G panel.
"""
import argparse
import concurrent.futures as cf
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

import _sumstats

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data" / "raw" / "gwas"
OUT = ROOT / "results" / "tables" / "02-mr"
CACHE = ROOT / "data" / "processed" / "mr" / "ld_cache"
FILES = _sumstats.FILES


def read_significant(key, p_thresh):
    """Genome-wide significant variants on the common schema, with F and R2."""
    d = _sumstats.read_rows(key, p_thresh)
    if d.empty:
        return d
    for c in ["pos", "eaf", "beta", "se", "p", "n"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["SNP", "beta", "se", "p", "eaf"])
    d = d[(d["eaf"] > 0.01) & (d["eaf"] < 0.99)]
    d = d[d["se"] > 0]
    d = d.sort_values("p").drop_duplicates("SNP", keep="first")
    d = d[d["effect_allele"].str.len() == 1]
    d = d[d["other_allele"].str.len() == 1]
    d["F"] = (d["beta"] / d["se"]) ** 2
    d["r2"] = d["F"] / (d["F"] + d["n"] - 2)
    return d.reset_index(drop=True)


def ld_pair(a, b, pop, tries=3):
    """r2 between two rsIDs from Ensembl REST (1000G phase 3, one population)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    tag = pop.split(":")[-1]
    cpath = CACHE / f"{tag}_{a}_{b}.json"
    if cpath.exists():
        try:
            return json.loads(cpath.read_text()).get("r2")
        except Exception:
            pass
    val = None
    for attempt in range(tries):
        try:
            url = (f"https://rest.ensembl.org/ld/human/pairwise/{a}/{b}"
                   f"?population_name={pop}")
            req = urllib.request.Request(url, headers={"Accept": "application/json",
                                                       "User-Agent": "research (module 02-mr)"})
            with urllib.request.urlopen(req, timeout=90) as r:
                recs = json.loads(r.read().decode())
            val = recs[0].get("r2") if recs else None
            break
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):      # not on the 1000G panel: no LD information
                break
            time.sleep(1.5 + 2 * attempt)
        except Exception:
            time.sleep(1.5 + 2 * attempt)
    cpath.write_text(json.dumps({"a": a, "b": b, "pop": pop, "r2": val}))
    return val


def prune_distance(d, kb):
    """Greedy distance pruning on the same chromosome (strongest p-value wins the window)."""
    kept, dropped = [], {}
    for chrom, g in d.sort_values("p").groupby("chr", sort=False):
        leads = []          # (pos, SNP) kept on this chromosome, sorted as we go
        for _, row in g.iterrows():
            hit = next((l for l in leads if abs(l[0] - row["pos"]) <= kb * 1000), None)
            if hit is None:
                leads.append((row["pos"], row["SNP"]))
                kept.append(row)
            else:
                dropped[row["SNP"]] = [hit[1], None]
    return pd.DataFrame(kept).sort_values("p").reset_index(drop=True), dropped


def clump_ld(d, r2_thresh, kb, pop, threads):
    """Greedy LD clumping: the strongest variant keeps the locus, its correlates are dropped."""
    leads, dropped, no_ld = [], {}, 0
    chrom = dict(zip(d["SNP"], d["chr"].astype(str)))
    pos = dict(zip(d["SNP"], d["pos"].astype("int64")))
    remaining = d["SNP"].tolist()
    n_req = 0
    with cf.ThreadPoolExecutor(max_workers=threads) as ex:
        while remaining:
            lead = remaining[0]
            near = [r for r in remaining[1:]
                    if chrom[r] == chrom[lead] and abs(pos[r] - pos[lead]) <= kb * 1000]
            futs = {ex.submit(ld_pair, lead, o, pop): o for o in near}
            n_req += len(futs)
            for f, other in futs.items():
                try:
                    r2 = f.result()
                except Exception:
                    r2 = None
                if r2 is None:
                    continue
                if r2 > r2_thresh and other in remaining:
                    remaining.remove(other)
                    dropped[other] = [lead, r2]
            if not near:
                no_ld += 1
            leads.append(lead)
            remaining.remove(lead)
    print(f"  {n_req} Ensembl LD pair queries")
    return leads, dropped, no_ld


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", default="eBMD")
    ap.add_argument("--p", type=float, default=5e-8)
    ap.add_argument("--clump", choices=["ld", "distance", "none"], default="ld")
    ap.add_argument("--r2", type=float, default=0.001)
    ap.add_argument("--kb", type=int, default=500)
    ap.add_argument("--pre-kb", type=int, default=100,
                    help="distance pre-pruning done before LD clumping, to bound LD queries")
    ap.add_argument("--pop", default="1000GENOMES:phase_3:CEU")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--keep-all", action="store_true",
                    help="also write the unclumped candidate list as instruments_<KEY>_all.csv "
                         "so that 02b_ld_clump.py can clump it with a local 1000G panel")
    a = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    summary = []
    for key in [k.strip() for k in a.datasets.split(",") if k.strip()]:
        path = RAW / FILES[key]
        if not path.exists():
            print(f"{key}: {path.name} not downloaded yet, skipped")
            continue
        print(f"\n=== {key}: {path.name}")
        t0 = time.time()
        d = read_significant(key, a.p)
        print(f"  {len(d)} variants with p < {a.p:g} (MAF 1-99%), {time.time() - t0:.0f}s")
        if d.empty:
            continue
        method = a.clump
        if a.clump == "ld":
            pre, pre_dropped = prune_distance(d, a.pre_kb) if a.pre_kb > 0 else (d, {})
            print(f"  distance pre-prune ({a.pre_kb} kb): {len(pre)} of {len(d)} kept")
            try:
                leads, dropped, no_ld = clump_ld(pre, a.r2, a.kb, a.pop, a.threads)
                dropped = {**pre_dropped, **dropped}
                print(f"  LD clumped ({a.pop}, r2>{a.r2}, {a.kb} kb): {len(leads)} independent "
                      f"(dropped {len(dropped)}, loci without LD info {no_ld})")
            except Exception as e:
                print(f"  Ensembl LD clumping failed ({e}); falling back to distance pruning")
                pre, dropped = prune_distance(d, a.kb)
                leads, no_ld = pre["SNP"].tolist(), 0
                method = "distance-fallback"
        elif a.clump == "distance":
            pre, dropped = prune_distance(d, a.kb)
            leads, no_ld = pre["SNP"].tolist(), 0
            print(f"  distance pruned ({a.kb} kb): {len(leads)} (dropped {len(dropped)})")
        else:
            leads, dropped, no_ld = d["SNP"].tolist(), {}, 0
            print(f"  no clumping: {len(leads)} variants (will over-count correlated signals)")
        inst = d[d["SNP"].isin(leads)].copy()
        cols = ["SNP", "chr", "pos", "effect_allele", "other_allele", "eaf", "beta", "se", "p", "n", "F", "r2"]
        inst = inst[cols].sort_values("p").reset_index(drop=True)
        if a.keep_all:
            # full candidate list (distance pre-pruned only) for local LD clumping
            cand = d[["SNP", "chr", "pos", "effect_allele", "other_allele", "eaf", "beta", "se", "p", "n", "F", "r2"]]
            cand.sort_values("p").reset_index(drop=True).to_csv(OUT / f"instruments_{key}_all.csv", index=False)
        inst.to_csv(OUT / f"instruments_{key}.csv", index=False)
        summary.append(dict(dataset=key, significant=len(d), instruments=len(inst), dropped=len(dropped),
                            no_ld_info=no_ld, clump_method=method, min_F=round(inst["F"].min(), 1),
                            mean_r2=round(inst["r2"].mean(), 6)))
        print(f"  -> {OUT.name}/instruments_{key}.csv  min F = {inst['F'].min():.1f}")
    if summary:
        s = pd.DataFrame(summary)
        s.to_csv(OUT / "instruments_summary.csv", index=False)
        print("\n" + s.to_string(index=False))


if __name__ == "__main__":
    main()

