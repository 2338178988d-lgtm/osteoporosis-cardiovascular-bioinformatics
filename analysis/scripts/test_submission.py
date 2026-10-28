"""Round-trip test: extract the package elsewhere and check it stands on its own.

A package that only works inside the project it was built from is not a package.
"""
import hashlib
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ZIP = next((ROOT / "results" / "submission").glob("*.zip"))
DEST = ROOT / ".cache" / "submission_test"

if DEST.exists():
    shutil.rmtree(DEST)
DEST.mkdir(parents=True)

print(f"extracting {ZIP.name} ({ZIP.stat().st_size / 1e6:.1f} MB) -> {DEST}")
with zipfile.ZipFile(ZIP) as zf:
    bad = zf.testzip()
    print(f"  zip integrity: {'OK' if bad is None else 'CORRUPT at ' + bad}")
    zf.extractall(DEST)

# 1. manifest verification, exactly as a reviewer would run it
man = DEST / "MANIFEST_SHA256.txt"
mismatch, missing = [], []
for line in man.read_text(encoding="utf-8").splitlines():
    if not line.strip():
        continue
    want, rel = line.split("  ", 1)
    f = DEST / rel
    if not f.exists():
        missing.append(rel)
        continue
    h = hashlib.sha256(f.read_bytes()).hexdigest()
    if h != want:
        mismatch.append(rel)
print(f"\nmanifest check: {len(man.read_text(encoding='utf-8').splitlines())} entries")
print(f"  missing   {len(missing)}  {missing[:3]}")
print(f"  mismatch  {len(mismatch)}  {mismatch[:3]}")

# 2. the manuscript verifier must run against the extracted tables
vm = DEST / "analysis" / "scripts" / "verify_manuscript.py"
r = subprocess.run([sys.executable, str(vm)], capture_output=True, text=True, cwd=DEST)
tail = (r.stdout or r.stderr).strip().splitlines()[-1] if (r.stdout or r.stderr) else ""
print(f"\nverify_manuscript.py from inside the package: exit {r.returncode}")
print(f"  {tail}")
if r.returncode != 0:
    print("  stderr tail:", (r.stderr or "")[-400:])

# 3. one figure must rebuild from inside the package
fig = DEST / "analysis" / "scripts" / "make_paper_fig1.py"
r2 = subprocess.run([sys.executable, str(fig)], capture_output=True, text=True, cwd=DEST)
lines = [x for x in (r2.stdout or "").splitlines() if "collision" in x or "saved" in x]
print(f"\nmake_paper_fig1.py from inside the package: exit {r2.returncode}")
for line in lines:
    print("  " + line)
if r2.returncode != 0:
    print("  stderr tail:", (r2.stderr or "")[-500:])

ok = not missing and not mismatch and r.returncode == 0 and r2.returncode == 0
print("\n" + ("PACKAGE IS SELF-CONTAINED" if ok else "PACKAGE HAS A PROBLEM"))

