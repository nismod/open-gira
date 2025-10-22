"""Check missing slices from a directory

Example usage:
    python check_missing_slices.py \
        results/direct_damages/planet-latest_filter-rail/hazard-landslide-arup/EAD_and_cost_per_trigger/

"""

import sys
from pathlib import Path
import pyarrow.parquet as pq
from tqdm import tqdm

dirn = sys.argv[1]
MAX_SLICE = 2025  # maybe add as a parameter? Current work is all with slice_count=2025
p = Path(dirn)

# Find and print the difference between expected and actual set of geoparquet slice files
exp = [f"slice-{n}.geoparquet" for n in range(MAX_SLICE)]
act = []

for fn in tqdm(list(p.glob("*.geoparquet"))):
    try:
        pq.read_schema(fn)
        act.append(fn.name)
    except:  # noqa: E722
        print("Error reading", fn)
        continue

nums = sorted(
    [
        int(fname.replace("slice-", "").replace(".geoparquet", ""))
        for fname in set(exp) - set(act)
    ]
)
if not nums:
    sys.exit()

prev = nums[0]
in_range = False
range_start = -1

k = []
start, stop = None, None
for num in nums:
    if stop is None:
        start, stop = num, num
    elif stop == num - 1:
        stop += 1
    else:
        if start == stop:
            k.append(str(start))
        else:
            k.append(f"{start}-{stop}")
        start, stop = num, num

if start == stop:
    k.append(str(start))
else:
    k.append(f"{start}-{stop}")

print(",".join(k))
