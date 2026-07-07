"""
Reader for STORM synthetic tracks (Bloemendaal et al.).

https://doi.org/10.1038/s41597-020-0381-2

STORM is distributed as headerless CSV, one file per basin per 1,000-year
sample, e.g. ``STORM_DATA_IBTRACS_EP_1000_YEARS_0.csv``.
"""

from os import PathLike
from typing import Iterable, Union

import geopandas as gpd
import numpy as np
import pandas as pd

from .. import ops
from ..core import TrackSet
from ..units import TEN_MINUTE_TO_ONE_MINUTE_WIND_FACTOR

#: Column names and dtypes of the raw STORM CSV files.
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

#: Basins are serialized as integers in the data: 0 -> "EP", 2 -> "NI", etc.
STORM_BASIN_IDS = ("EP", "NA", "NI", "SI", "SP", "WP")

#: Temporal frequency of STORM synthetic track points.
STORM_FREQUENCY = "3h"

#: Simulated years per STORM sample file. Basins are simulated concurrently,
#: so files for different basins of the same sample share one 1,000 year span.
STORM_YEARS_PER_SAMPLE = 1_000.0


def read_storm(
    paths: Iterable[Union[str, PathLike]],
    sample: int,
    source: str = "STORM-constant",
) -> TrackSet:
    """
    Read STORM CSV files (one or more basins) for a single sample.

    Wind speeds are rescaled from 10-minute to 1-minute sustained. The
    ``year`` column is offset by ``1000 * sample`` so that concatenating
    samples yields unique (year, tc_number) pairs, matching open-gira's
    existing convention.
    """

    per_basin = []
    for path in paths:
        df = pd.read_csv(path, names=STORM_CSV_SCHEMA.keys(), dtype=STORM_CSV_SCHEMA)

        df["sample"] = int(sample)
        df["year"] = df["year"].astype(int) + 1_000 * int(sample)
        df["lon"] = ops.wrap_longitude(df["lon"])

        # lookup string basin code from integer representation
        df["basin_id"] = np.array(STORM_BASIN_IDS)[df["basin_id"]]

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
        df = ops.synthetic_time_index(df, freq=STORM_FREQUENCY)

        per_basin.append(df)

    df = ops.concat_tracks(per_basin)

    # rescale winds from 10-minutely to 1-minutely
    df["max_wind_speed_ms"] /= TEN_MINUTE_TO_ONE_MINUTE_WIND_FACTOR

    df = gpd.GeoDataFrame(
        data=df.drop(columns=["lat", "lon"]),
        geometry=gpd.points_from_xy(df["lon"], df["lat"], crs=4326),
    )

    return TrackSet(
        data=df,
        source=source,
        years=STORM_YEARS_PER_SAMPLE,
        synthetic_time=True,
    )
