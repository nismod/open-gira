"""Check missing slices from a directory

Example usage:
    python check_missing.py results/direct_damages/planet-latest_filter-rail/hazard-landslide-arup/EAD_and_cost_per_trigger/

"""

import sys
from pathlib import Path

dirn = sys.argv[1]
MAX_SLICE = 2025  # maybe add as a parameter? Current work is all with slice_count=2025
p = Path(dirn)

# Find and print the difference between expected and actual set of geoparquet slice files
exp = [f"slice-{n}.geoparquet" for n in range(MAX_SLICE)]
act = [fn.name for fn in p.glob("*.geoparquet")]
for fn in sorted(list(set(exp) - set(act))):
    print(fn)
