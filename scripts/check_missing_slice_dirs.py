"""Check missing slices from a directory

Example usage:
    python check_missing_slice_dirs.py results/slices/planet-latest_nbs/

"""

import sys
from pathlib import Path
import pyarrow.parquet as pq
from tqdm import tqdm

dirn = sys.argv[1]
fname = sys.argv[2]
MAX_SLICE = 2025  # maybe add as a parameter? Current work is all with slice_count=2025
p = Path(dirn)

# Find and print the difference between expected and actual set of geoparquet slice files
exp = [f"slice-{n}" for n in range(MAX_SLICE)]
act = []

for fn in tqdm(list(p.glob(f"slice-*/{fname}"))):
    try:
        pq.read_schema(fn)
        act.append(fn.parent.name)
    except:
        print("Error reading", fn)
        continue

nums = sorted([int(fname.replace("slice-", "")) for fname in set(exp) - set(act)])
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
