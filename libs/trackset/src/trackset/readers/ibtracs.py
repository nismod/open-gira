"""
Reader for IBTrACS observed best-track data (Knapp et al.).

https://www.ncei.noaa.gov/products/international-best-track-archive

Column documentation:
https://www.ncei.noaa.gov/sites/default/files/2021-07/IBTrACS_v04_column_documentation.pdf
"""

from os import PathLike
from typing import Union

import geopandas as gpd
import pandas as pd

from .. import ops
from ..core import TrackSet
from ..units import (
    IBTRACS_AGENCY_1MIN_WIND_FACTOR,
    KM_PER_NAUTICAL_MILE,
    MS_PER_KNOT,
)

_FIXED_COLUMNS = {
    "ISO_TIME": str,
    "LON": float,
    "LAT": float,
    "SID": str,
    "NAME": str,
    "SEASON": int,  # year
    "NUMBER": int,  # cardinal number of storm in year
    "BASIN": str,
    "DIST2LAND": float,  # km to nearest continent or island larger than 1400km^2
}


def read_ibtracs(path: Union[str, PathLike]) -> TrackSet:
    """
    Read an IBTrACS CSV file.

    Wind speeds, radii and pressures reported by multiple agencies are
    normalised (winds to a 1-minute averaging period) and averaged across
    agencies. Points missing any of the wind, pressure or radius variables
    required for wind field estimation are dropped.
    """

    # read the header line to assemble the lists of per-agency wind and
    # pressure columns present in this version of the file
    with open(path, "r") as fp:
        header: list[str] = fp.readline().strip().split(",")

    wind_cols = [c for c in header if c.endswith("_WIND")]  # knots
    rmw_cols = [c for c in header if c.endswith("_RMW")]  # nautical miles
    pressure_cols = [c for c in header if c.endswith("_PRES")]  # mb / hPa

    df = pd.read_csv(
        path,
        skiprows=[1],  # drop units definition row
        usecols=list(_FIXED_COLUMNS) + wind_cols + rmw_cols + pressure_cols,
        dtype=_FIXED_COLUMNS
        | {col: float for col in wind_cols + rmw_cols + pressure_cols},
        header=0,
        keep_default_na=False,  # otherwise 'NA' (North Atlantic) is read as NaN!
        na_values=["", " "],  # missing data value for IBTrACS CSV is a space
    )
    df = df.rename(
        columns={
            "ISO_TIME": "time_utc",
            "LON": "lon",
            "LAT": "lat",
            "SID": "track_id",
            "NAME": "name",
            "SEASON": "year",
            "NUMBER": "tc_number",
            "BASIN": "basin_id",
            "DIST2LAND": "distance_to_land_km",
        }
    )

    # normalise each agency's winds to a 1-minute averaging period
    for agency, (scale, shift) in IBTRACS_AGENCY_1MIN_WIND_FACTOR.items():
        col = f"{agency}_WIND"
        if col in df.columns:
            df[col] = (df[col] - shift) / scale

    # average reports across agencies
    df["max_wind_speed_ms"] = (
        df.loc[:, wind_cols].mean(axis="columns", skipna=True) * MS_PER_KNOT
    )
    # there are some tens of negative valued wind observations
    df.loc[df["max_wind_speed_ms"] < 0, "max_wind_speed_ms"] = 0.0

    df["radius_to_max_winds_km"] = (
        df.loc[:, rmw_cols].mean(axis="columns", skipna=True) * KM_PER_NAUTICAL_MILE
    )

    df["min_pressure_hpa"] = df.loc[:, pressure_cols].mean(axis="columns", skipna=True)
    # there are some hundred or so of implausibly low pressure observations
    df.loc[df["min_pressure_hpa"] < 800, "min_pressure_hpa"] = 1000.0

    df["category"] = ops.saffir_simpson_category(df["max_wind_speed_ms"])

    # keep only points with the variables necessary for wind field estimation
    df = df[
        df["max_wind_speed_ms"].notna()
        & df["min_pressure_hpa"].notna()
        & df["radius_to_max_winds_km"].notna()
    ].copy()
    df["category"] = df["category"].astype(int)

    df["timestep"] = df.groupby("track_id", sort=False).cumcount()

    df = df.set_index(pd.to_datetime(df["time_utc"])).drop(columns=["time_utc"])
    df.index.name = None
    df["year"] = df.index.year
    df["month"] = df.index.month

    # boolean: is storm over land?
    df["landfall"] = df["distance_to_land_km"] == 0

    df = df.drop(columns=wind_cols + rmw_cols + pressure_cols)

    df["lon"] = ops.wrap_longitude(df["lon"])
    df = gpd.GeoDataFrame(
        data=df.drop(columns=["lat", "lon"]),
        geometry=gpd.points_from_xy(df["lon"], df["lat"], crs=4326),
    )

    # observed record: the set spans the years present in the data
    years = float(df["year"].max() - df["year"].min() + 1)

    return TrackSet(
        data=df,
        source="IBTrACS",
        years=years,
        synthetic_time=False,
    )
