"""
Concatenate CSV files for all IRIS synthetic tracks of a given model
into a single geoparquet file.

Output format shares columns with IBTrACS.
"""

from glob import glob
import os
import re
from typing import List

import numpy as np
import geopandas as gpd
import pandas as pd
from tqdm import tqdm

from open_gira.utils import natural_sort


# column names and dtypes for IRIS synthetic tropical cyclone tracks
IRIS_CSV_SCHEMA = {
    "year": int,
    "tc": int,
    "month": int,
    "timestep": int,
    "lon": float,
    "lat": float,
    "vmax": float,
    "pmin": float,
    "rmw": float,
    "r18": float
}
# basins are serialized as integers in data, 0 -> "SI", 2 -> "WP" etc.
IRIS_BASIN_IDS = ("SI", "SP", "WP", "EP", "NA", "NI")
# temporal frequency of IRIS
IRIS_FREQUENCY = "3H"


if __name__ == "__main__":

    csv_dir = snakemake.input.csv_dir
    parquet_path = snakemake.output.parquet

    data = []
    for path in tqdm(natural_sort(glob(f"{csv_dir}/*.txt"))):

        df = pd.read_csv(path, header=1, sep=r"\s+", names=IRIS_CSV_SCHEMA.keys(), dtype=IRIS_CSV_SCHEMA)
        df = df.rename(
            columns={
                "tc": "tc_number",
                "vmax": "max_wind_speed_ms",
                "rmw": "radius_to_max_winds_km",
                "pmin": "min_pressure_hpa",
            }
        )
        df = df.drop(columns=["r18"])

        # example paths containing sample number and basin id:
        # IRIS_PRESENT_bas0_1000Y_n1.txt
        # IRIS_2050-SSP1_bas0_1000Y_n6.txt
        sample, = re.search(r"1000Y_n([\d]).txt", os.path.basename(path)).groups()
        basin_id, = re.search(r"_bas([\d])_1000Y_", os.path.basename(path)).groups()

        df["sample"] = int(sample)
        df["basin_id"] = int(basin_id)

        # lookup string basin code from integer representation
        df.basin_id = np.array(IRIS_BASIN_IDS)[df.basin_id]

        # this gives us 10,000 years of data rather than 10 * 10,000 years
        # it is necessary when calculating the expected annual exposure/disruption to know the total timespan
        df["year"] = df["year"].astype(int) + 1000 * int(sample)

        # change geometry from 0-360 to -180-180
        df.lon = np.where(df.lon > 180, df.lon - 360, df.lon)

        df["track_id"] = (
            df["basin_id"] + "_"
            + df["sample"].astype(str) + "_"
            + df["year"].astype(int).astype(str) + "_"
            + df["tc_number"].astype(int).astype(str)
        )

        # we'll want to interpolate and then measure the speed of tracks later,
        # this is easiest when we have some temporal index (as in IBTrACS)
        # so make up an artificial one here based on the IRIS reporting frequency

        track_datetimes: List[np.ndarray] = []
        track_lengths: np.ndarray = df.track_id.apply(hash).value_counts(sort=False).values
        for length in track_lengths:
            track_datetimes.append(pd.date_range(start="2000-01-01", periods=length, freq=IRIS_FREQUENCY).values)

        df = df.set_index(np.concatenate(track_datetimes))

        data.append(df)

    df = pd.concat(data)

    # construct geometry from lat and long
    df = gpd.GeoDataFrame(
        data=df,
        geometry=gpd.points_from_xy(df["lon"], df["lat"], crs=4326)
    )
    df = df.drop(columns=["lon", "lat"])

    os.makedirs(os.path.dirname(parquet_path), exist_ok=True)
    df.to_parquet(parquet_path)
