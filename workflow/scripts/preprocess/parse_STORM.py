"""
Concatenate CSV files for all STORM synthetic tracks of a given model
into a single geoparquet file.

Output format shares columns with IBTrACS.
"""

from glob import glob
import logging
import os
import re
from typing import List

import numpy as np
import geopandas as gpd
import pandas as pd
from tqdm import tqdm

from open_gira.utils import natural_sort


# column names and dtypes for synthetic tropical cyclone tracks
STORM_CSV_SCHEMA = {
    "year": int,
    "month": int,
    "tc_number": int,
    "timestep": int,
    "basin_id": int,
    "lat": float,
    "lon": float,
    "min_pressure_hpa": float,
    "max_wind_speed_ms": float,
    "radius_to_max_winds_km": float,
    "category": int,
    "landfall": bool,  # actually {0|1} in the file
    "distance_to_land_km": float,
}
# basins are serialized as integers in data, 0 -> "EP", 2 -> "NI" etc.
STORM_BASIN_IDS = ("EP", "NA", "NI", "SI", "SP", "WP")
# divide by this factor to 'convert' STORM's 10-minutely sustained winds to
# 1-minutely sustained wind speeds, noting the vagueries of this process as
# explained here: https://library.wmo.int/doc_num.php?explnum_id=290
STORM_1MIN_WIND_FACTOR = 0.88
# temporal frequency of STORM synthetic tracks
STORM_FREQUENCY = "3H"


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    csv_dir = snakemake.input.csv_dir
    parquet_path = snakemake.output.parquet
    sample = snakemake.wildcards.SAMPLE

    data = []
    # loop over basins for given sample, processing tracks into a common format
    for path in natural_sort(glob(f"{csv_dir}/*_1000_YEARS_{sample}*.csv")):

        logging.info(path)

        df = pd.read_csv(path, names=STORM_CSV_SCHEMA.keys(), dtype=STORM_CSV_SCHEMA)

        df["sample"] = int(sample)

        # this gives us 10,000 years of data rather than 10 * 10,000 years
        # it is necessary when calculating the expected annual exposure/disruption to know the total timespan
        df["year"] = df["year"].astype(int) + 1000 * int(sample)

        # change geometry from 0-360 to -180-180
        df.lon = np.where(df.lon > 180, df.lon - 360, df.lon)

        # lookup string basin code from integer representation
        df.basin_id = np.array(STORM_BASIN_IDS)[df.basin_id]

        # different track_id format for STORM vs. IBTrACS, ensures no collisions
        df["track_id"] = (
            df["basin_id"] + "_"
            + df["sample"].astype(str) + "_"
            + df["year"].astype(int).astype(str) + "_"
            + df["tc_number"].astype(int).astype(str)
        )

        # STORM contains duplicate track points (duplicate timesteps for a given year-tc_number combination)
        df["point_id"] = df.apply(lambda row: f"{row.track_id}_{row.timestep}", axis=1)
        n_rows_raw = len(df)
        df = df.drop_duplicates(subset="point_id").drop(columns=["point_id"])
        logging.info(f"Collated {n_rows_raw} track points, dropped {n_rows_raw - len(df)} as duplicates")

        # we'll want to interpolate and then measure the speed of tracks later,
        # this is easiest when we have some temporal index (as in IBTrACS)
        # so make up an artificial one here based on the STORM reporting frequency
        track_datetimes: List[np.ndarray] = []
        track_lengths: np.ndarray = df.track_id.apply(hash).value_counts(sort=False).values
        for length in track_lengths:
            track_datetimes.append(pd.date_range(start="2000-01-01", periods=length, freq=STORM_FREQUENCY).values)

        df = df.set_index(np.concatenate(track_datetimes))

        # reorder columns
        df = df.loc[:, list(STORM_CSV_SCHEMA.keys()) + ["track_id", "sample"]]

        data.append(df)

    df = pd.concat(data)

    # rescale winds to 1-minutely
    df.max_wind_speed_ms /= STORM_1MIN_WIND_FACTOR

    # construct geometry from lat and long
    df = gpd.GeoDataFrame(
        data=df,
        geometry=gpd.points_from_xy(df["lon"], df["lat"], crs=4326)
    )
    df = df.drop(columns=["lon", "lat"])

    os.makedirs(os.path.dirname(parquet_path), exist_ok=True)
    df.to_parquet(parquet_path)
