"""
Concatenate CSV files for all STORM synthetic tracks of a given model
into a single geoparquet file.

Output format shares columns with IBTrACS.
"""

import os
import re
from glob import glob

import numpy as np
import geopandas as gpd
import pandas as pd

from open_gira.io import STORM_BASIN_IDS
from open_gira.io import STORM_CSV_SCHEMA as schema
from open_gira.utils import natural_sort


# divide by this factor to 'convert' STORM's 10-minutely sustained winds to
# 1-minutely sustained wind speeds, noting the vagueries of this process as
# explained here: https://library.wmo.int/doc_num.php?explnum_id=290
STORM_1MIN_WIND_FACTOR = 0.88


if __name__ == "__main__":

    csv_dir = snakemake.input.csv_dir
    parquet_path = snakemake.output.parquet

    data = []
    for path in natural_sort(glob(f"{csv_dir}/*.csv")):

        df = pd.read_csv(path, names=schema.keys(), dtype=schema)

        # example paths containing sample number:
        # STORM_DATA_HadGEM3-GC31-HM_WP_1000_YEARS_9_IBTRACSDELTA.csv
        # STORM_DATA_IBTRACS_EP_1000_YEARS_0.csv
        sample, = re.search(r"1000_YEARS_([\d])", os.path.basename(path)).groups()

        df["sample"] = int(sample)

        # different track_id format for STORM vs. IBTrACS, ensures no collisions
        df["track_id"] = (
            df["sample"].astype(str)
            + "_"
            + df["year"].astype(int).astype(str)
            + "_"
            + df["tc_number"].astype(int).astype(str)
        )

        # change geometry from 0-360 to -180-180
        df.lon = np.where(df.lon > 180, df.lon - 360, df.lon)

        # lookup string basin code from integer representation
        df.basin_id = np.array(STORM_BASIN_IDS)[df.basin_id]

        # reorder columns
        df = df.loc[:, list(schema.keys()) + ["track_id", "sample"]]

        data.append(df)

    df = pd.concat(data).reset_index(drop=True)

    # rescale winds to 1-minutely
    df.max_wind_speed_ms /= STORM_1MIN_WIND_FACTOR

    # construct geometry from lat and long
    df = gpd.GeoDataFrame(
        data=df,
        geometry=gpd.points_from_xy(df["lon"], df["lat"], crs=4326)
    )
    df = df.drop(columns=["lon", "lat"])

    df.to_parquet(parquet_path)
