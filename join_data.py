import sys
from os.path import join, dirname
import geopandas as gpd


def append_data(base, slice_files):
    slice_files.pop()
    if len(slice_files) == 0:
        return
    base = base.append(gpd.read_parquet(slice_files[-1]))
    append_data(base, slice_files)


slice_files = sys.argv[1:]
slice_files = slice_files[::-1]

base = gpd.read_parquet(slice_files[-1])
append_data(base, slice_files)
base.to_parquet(
    join(
        "outputs",
        "tanzania-latest.splits.geoparquet"
    )
)
