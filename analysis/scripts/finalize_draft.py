"""Final consistency pass over the manuscript draft."""
import re
from pathlib import Path

t = Path("results/reports/璁烘枃鑽夌_v1.md").read_text(encoding="utf-8")

# 1. every figure referenced exists
figs = set(re.findall(r"figures/paper/(fig\d_[a-z_]+\.png)", t))
have = {p.name for p in Path("figures/paper").glob("*.png")}
print("figures referenced:", sorted(figs))
print("missing:", sorted(figs - have) or "none")
print("present but unreferenced:", sorted(have - figs) or "none")

# 2. section structure
print("\nsections:")
for line in t.splitlines():
    if re.match(r"^#{1,3} ", line):
        print("  " + line)

# 3. no leftover placeholders / stray escapes
for bad in ["\\\\n", "TODO", "TBD", "XXX", "寰呰ˉ", "???"]:
    if bad in t:
        print(f"\n!! found {bad!r}")

# 4. character count by section
parts = re.split(r"^## ", t, flags=re.M)
print("\nlength by top-level section:")
for p in parts[1:]:
    name = p.splitlines()[0]
    print(f"  {name:24s} {len(p):6d} chars")
print(f"  {'TOTAL':24s} {len(t):6d} chars")

