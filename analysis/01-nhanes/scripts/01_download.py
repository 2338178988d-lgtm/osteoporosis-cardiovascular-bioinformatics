"""Step 1 鈥?download NHANES components needed for the osteoporosis x CVD analysis.

Files are fetched from the current CDC URL scheme
(https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/<year>/DataFiles/<FILE>.XPT).
That scheme answers 200 with text/html for files that do not exist in a cycle, so
availability is decided on Content-Type + magic bytes, never on status alone.

Outputs
  data/raw/nhanes/<cycle>/<FILE>.XPT      raw, immutable
  data/raw/nhanes/manifest.csv            url, bytes, sha256, rows per file
"""
import csv
import hashlib
import io
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

RAW = Path("data/raw/nhanes")
CYCLES = [("2005-2006", 2005, "D"), ("2007-2008", 2007, "E"), ("2009-2010", 2009, "F"),
          ("2011-2012", 2011, "G"), ("2013-2014", 2013, "H"), ("2015-2016", 2015, "I"),
          ("2017-2018", 2017, "J")]
# component name -> why it is needed
WANTED = {
    "DEMO": "age, sex, race/ethnicity, income, education, survey design",
    "BMX": "BMI, waist",
    "DXXFEM": "femoral neck / total femur / trochanter BMD (DXA), split-file cycles",
    "DXXSPN": "lumbar spine BMD (DXA), split-file cycles",
    "DXX": "combined DXA file (femur+spine) used in 2011-2012 and 2015-2016",
    "OSQ": "self-reported osteoporosis diagnosis",
    "MCQ": "self-reported CHF, CHD, angina, MI, stroke",
    "BPQ": "hypertension diagnosis and treatment",
    "BPX": "measured blood pressure (auscultatory / oscillometric)",
    "DIQ": "diabetes diagnosis",
    "SMQ": "smoking",
    "ALQ": "alcohol",
    "TCHOL": "total cholesterol, HDL (2005-2014)",
    "HDL": "HDL (2015-2018)",
    "TRIGLY": "triglycerides (2015-2018)",
    "BIOPRO": "creatinine, albumin, calcium, phosphorus, ALP",
    "GHB": "HbA1c",
    "RHQ": "reproductive history (menopause)",
}
UA = {"User-Agent": "Mozilla/5.0 (wisp-science research data acquisition)"}
MORT_BASE = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality"


def url_for(year, fname):
    return f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{year}/DataFiles/{fname}.XPT"


def fetch(url, timeout=180):
    """Return bytes for a genuine XPORT file, else None (soft-404 HTML is rejected)."""
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
                ctype = (r.headers.get("Content-Type") or "").lower()
                data = r.read()
            if "html" in ctype or len(data) < 2000:
                return None
            if b"HEADER RECORD" not in data[:400]:
                return None
            return data
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            time.sleep(2 * (attempt + 1))
        except Exception:
            time.sleep(2 * (attempt + 1))
    return None


def main():
    rows = []
    for cycle, year, suffix in CYCLES:
        outdir = RAW / cycle
        outdir.mkdir(parents=True, exist_ok=True)
        for comp, why in WANTED.items():
            fname = f"{comp}_{suffix}"
            dest = outdir / f"{fname}.XPT"
            if dest.exists() and dest.stat().st_size > 2000:
                data = dest.read_bytes()
                status = "cached"
            else:
                data = fetch(url_for(year, fname))
                status = "downloaded"
                if data is None:
                    print(f"  - {cycle} {fname:10s} not available in this cycle")
                    continue
                dest.write_bytes(data)
                time.sleep(0.5)
            try:
                n = len(pd.read_sas(io.BytesIO(data), format="xport"))
            except Exception:
                n = -1
            rows.append({
                "cycle": cycle, "component": comp, "file": f"{fname}.XPT",
                "url": url_for(year, fname), "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(), "n_rows": n,
                "status": status, "purpose": why,
            })
            print(f"  + {cycle} {fname:10s} rows={n:>6} bytes={len(data):>9} [{status}]")

    # linked mortality (fixed-width public-use files, 2019 follow-up)
    mortdir = RAW / "mortality"
    mortdir.mkdir(parents=True, exist_ok=True)
    for cycle in [c[0] for c in CYCLES]:
        fname = f"NHANES_{cycle.replace('-', '_')}_MORT_2019_PUBLIC.dat"
        dest = mortdir / fname
        if dest.exists() and dest.stat().st_size > 10000:
            data = dest.read_bytes()
            status = "cached"
        else:
            try:
                with urllib.request.urlopen(urllib.request.Request(f"{MORT_BASE}/{fname}", headers=UA), timeout=180) as r:
                    data = r.read()
                dest.write_bytes(data)
                status = "downloaded"
            except Exception as e:
                print(f"  - mortality {fname}: {e}")
                continue
        rows.append({"cycle": cycle, "component": "MORTALITY", "file": fname,
                     "url": f"{MORT_BASE}/{fname}", "bytes": len(data),
                     "sha256": hashlib.sha256(data).hexdigest(),
                     "n_rows": len(data.splitlines()), "status": status,
                     "purpose": "linked mortality through 2019-12-31 (all-cause, CVD)"})
        print(f"  + mortality {fname} lines={len(data.splitlines()):>6}")

    man = RAW / "manifest.csv"
    with open(man, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["cycle", "component", "file", "url", "bytes", "sha256", "n_rows", "status", "purpose"])
        w.writeheader()
        w.writerows(rows)
    print(f"\n{len(rows)} files -> {man}")


if __name__ == "__main__":
    sys.exit(main())

