"""
Concatenate CSV files for all STORM synthetic tracks of a given model
into a single geoparquet file.

Output format shares columns with IBTrACS.
"""

import os
from glob import glob

import numpy as np
import geopandas as gpd
import pandas as pd

from open_gira.io import STORM_BASIN_IDS
from open_gira.io import STORM_CSV_SCHEMA as schema
from open_gira.utils import natural_sort

if __name__ == "__main__":

    csv_dir = snakemake.input.csv_dir
    parquet_path = snakemake.output.parquet

    data = []
    for path in natural_sort(glob(f"{csv_dir}/*.csv")):
        sample = int(path.split(".csv")[0].split("_")[-1])
        df = pd.read_csv(path, names=schema.keys(), dtype=schema)

        df["sample"] = sample

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

    # construct geometry from lat and long
    df = gpd.GeoDataFrame(
        data=df,
        geometry=gpd.points_from_xy(df["lon"], df["lat"], crs=4326)
    )
    df = df.drop(columns=["lon", "lat"])

    df.to_parquet(parquet_path)
