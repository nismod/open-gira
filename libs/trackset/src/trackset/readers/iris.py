"""
Reader for IRIS synthetic tracks (Sparks & Toumi).

https://doi.org/10.1038/s41597-024-03250-y

IRIS is distributed as whitespace-separated text, one file per basin per
1,000-year sample, e.g. ``IRIS_WP_1000Y_n3.txt``, with two header lines.
"""

import os
import re
from os import PathLike
from typing import Iterable, Union

import geopandas as gpd
import pandas as pd

from .. import ops
from ..core import TrackSet

#: Column names and dtypes of the raw IRIS text files.
IRIS_CSV_SCHEMA = {
    "tcid": str,
    "year": int,
    "tc": int,
    "month": int,
    "timestep": int,
    "lon": float,
    "lat": float,
    "vmax": float,
    "pmin": float,
    "rmw": float,
    "r18": float,
}

#: Temporal frequency of IRIS synthetic track points.
IRIS_FREQUENCY = "3h"

#: Simulated years per IRIS sample file.
IRIS_YEARS_PER_SAMPLE = 1_000.0


def _sample_and_basin_from_filename(path: Union[str, PathLike]) -> tuple[int, str]:
    """
    Parse sample number and basin from names like ``IRIS_WP_1000Y_n3.txt``.
    """

    filename = os.path.basename(str(path))
    (sample,) = re.search(r"1000Y_n(\d+)", filename).groups()
    (basin_id,) = re.search(r"IRIS_(\w\w)_1000Y", filename).groups()
    return int(sample), basin_id


def read_iris(
    paths: Iterable[Union[str, PathLike]],
    source: str = "IRIS-PRESENT",
) -> TrackSet:
    """
    Read IRIS text files (one or more basins) for a single sample.

    Sample number and basin are parsed from the filenames. The ``year``
    column is offset by ``1000 * sample``, matching open-gira's existing
    convention (see :func:`trackset.readers.storm.read_storm`).
    """

    per_basin = []
    sample = None
    for path in paths:
        df = pd.read_csv(
            path,
            header=1,
            sep=r"\s+",
            names=IRIS_CSV_SCHEMA.keys(),
            dtype=IRIS_CSV_SCHEMA,
        )
        df = df.rename(
            columns={
                "tc": "tc_number",
                "vmax": "max_wind_speed_ms",
                "rmw": "radius_to_max_winds_km",
                "pmin": "min_pressure_hpa",
            }
        )
        df = df.drop(columns=["tcid", "r18"])

        sample, basin_id = _sample_and_basin_from_filename(path)
        df["sample"] = sample
        df["basin_id"] = basin_id
        df["year"] = df["year"].astype(int) + 1_000 * sample
        df["lon"] = ops.wrap_longitude(df["lon"])

        df["track_id"] = (
            df["basin_id"]
            + "_"
            + df["sample"].astype(str)
            + "_"
            + df["year"].astype(str)
            + "_"
            + df["tc_number"].astype(str)
        )

        df = ops.drop_duplicate_points(df)
        df = ops.synthetic_time_index(df, freq=IRIS_FREQUENCY)

        per_basin.append(df)

    df = ops.concat_tracks(per_basin)

    df = gpd.GeoDataFrame(
        data=df.drop(columns=["lat", "lon"]),
        geometry=gpd.points_from_xy(df["lon"], df["lat"], crs=4326),
    )

    return TrackSet(
        data=df,
        source=source,
        years=IRIS_YEARS_PER_SAMPLE,
        synthetic_time=True,
    )
