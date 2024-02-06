"""
Parse and clean IBTrACS historic storm CSV records into geoparquet

The output schema is a superset of the STORM synthetic track data
"""

import logging
import os
from typing import Union

import geopandas as gpd
import numpy as np
import pandas as pd

from parse_STORM import STORM_CSV_SCHEMA


"""
Wind scale and shift values taken from CLIMADA: climada.hazard.tc_tracks

We use these to 'convert' between the agencies' reporting periods and a
1-minutely period.

From Table 1 in: Knapp, K.R., Kruk, M.C. (2010): Quantifying Interagency
Differences in Tropical Cyclone Best-Track Wind Speed Estimates.
Monthly Weather Review 138(4): 1459–1473.
https://library.wmo.int/index.php?lvl=notice_display&id=135
"""
IBTRACS_AGENCY_1MIN_WIND_FACTOR = {
    "USA": [1.0, 0.0],
    "TOKYO": [0.60, 23.3],
    "NEWDELHI": [1.0, 0.0],
    "REUNION": [0.88, 0.0],
    "BOM": [0.88, 0.0],
    "NADI": [0.88, 0.0],
    "WELLINGTON": [0.88, 0.0],
    "CMA": [0.871, 0.0],
    "HKO": [0.9, 0.0],
    "DS824": [1.0, 0.0],
    "TD9636": [1.0, 0.0],
    "TD9635": [1.0, 0.0],
    "NEUMANN": [0.88, 0.0],
    "MLC": [1.0, 0.0],
}


def saffir_simpson_classifier(wind_speed_ms: float) -> Union[int, float]:
    """
    Identify the Saffir-Simpson storm category given a wind speed in m/s.

    The input wind speed should be a 1-minute sustained measurement
    """

    if not isinstance(wind_speed_ms, (float, int)):
        raise ValueError(f"{wind_speed_ms=} must be float or integer")

    if wind_speed_ms < 0:
        raise ValueError(f"{wind_speed_ms=} should be positive-valued")
    elif np.isnan(wind_speed_ms):
        return np.nan
    elif wind_speed_ms < 18:
        return -1  # Tropical Disturbance
    elif wind_speed_ms < 33:
        return 0  # Tropical Storm
    elif wind_speed_ms < 43:
        return 1
    elif wind_speed_ms < 50:
        return 2
    elif wind_speed_ms < 58:
        return 3
    elif wind_speed_ms < 70:
        return 4
    elif wind_speed_ms >= 70:
        return 5


if __name__ == "__main__":

    ibtracs_csv_path = snakemake.input.ibtracs_csv
    ibtracs_parquet_path = snakemake.output.ibtracs_parquet

    # read the header line to assemble a list of wind and pressure columns to read
    # https://www.ncei.noaa.gov/sites/default/files/2021-07/IBTrACS_v04_column_documentation.pdf
    with open(ibtracs_csv_path, "r") as fp:
        columns: list[str] = fp.readline().split(",")

    WIND_COLS = list(filter(lambda s: s.endswith("_WIND"), columns))  # maximum wind speeds in knots
    RMW_COLS = list(filter(lambda s: s.endswith("_RMW"), columns))  # radius from eye to maximum winds in nautical miles
    PRESSURE_COLS = list(filter(lambda s: s.endswith("_PRES"), columns))  # eye / minimum pressures in mb / hPa
    OTHER_COLS = {
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

    df = pd.read_csv(
        ibtracs_csv_path,
        skiprows=[1],  # drop units definition row
        usecols=list(OTHER_COLS.keys()) + WIND_COLS + RMW_COLS + PRESSURE_COLS,
        dtype=OTHER_COLS | {col: float for col in WIND_COLS + RMW_COLS + PRESSURE_COLS},
        header=0,
        keep_default_na=False,  # otherwise 'NA' (North Atlantic basin) is treated as a NaN value!
        na_values=['', ' ']  # missing data value for IBTrACS CSV is a single space
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

    for agency, (scale, shift) in IBTRACS_AGENCY_1MIN_WIND_FACTOR.items():
        df[f'{agency}_WIND'] -= shift
        df[f'{agency}_WIND'] /= scale

    # average wind columns from various agencies
    METERS_PER_SECOND_PER_KNOT = 0.5144
    df["max_wind_speed_ms"] = df.loc[:, WIND_COLS].apply(np.nanmean, axis="columns") * METERS_PER_SECOND_PER_KNOT
    # there are some tens of negative valued wind observations
    df.loc[df.loc[:, "max_wind_speed_ms"] < 0, "max_wind_speed_ms"] = 0

    # label with Saffir-Simpson storm category
    df["category"] = df["max_wind_speed_ms"].apply(saffir_simpson_classifier)

    # average radius values from agencies
    KM_PER_NAUTICAL_MILE = 1.852
    df["radius_to_max_winds_km"] = df.loc[:, RMW_COLS].apply(np.nanmean, axis="columns") * KM_PER_NAUTICAL_MILE

    # again, average values from agencies
    df["min_pressure_hpa"] = df.loc[:, PRESSURE_COLS].apply(np.nanmean, axis="columns")
    # there are some hundred or so of zero valued pressure observations
    df.loc[df.loc[:, "min_pressure_hpa"] < 800, "min_pressure_hpa"] = 1000

    # filter dataset by presence of variables necessary for computing wind field
    df = df[~df["max_wind_speed_ms"].isna() & ~df["min_pressure_hpa"].isna() & ~df["radius_to_max_winds_km"].isna()]

    # rest of processing is to match STORM synthetic track output

    for track_id in set(df.track_id):
        mask = (df.track_id == track_id)
        df.loc[mask, "timestep"] = range(len(df[mask]))
    df["timestep"] = df["timestep"].astype(int)

    df = df.set_index(pd.to_datetime(df["time_utc"]))
    df["year"] = df.index.year
    df["month"] = df.index.month

    # boolean: is storm over land?
    df["landfall"] = df["distance_to_land_km"] == 0

    # the synthetic tracksets have a sample attribute, match for common interface
    df["sample"] = 0

    df = df.drop(columns=WIND_COLS + RMW_COLS + PRESSURE_COLS + ["time_utc"])

    # reorder columns to match STORM dataset (except track_id and name, which are additional)
    df = df.loc[:, list(STORM_CSV_SCHEMA.keys()) + ["track_id", "name"]]

    # change geometry from 0-360 to -180-180
    df.lon = np.where(df.lon > 180, df.lon - 360, df.lon)

    # construct geometry from lat and long
    df = gpd.GeoDataFrame(
        data=df,
        geometry=gpd.points_from_xy(df["lon"], df["lat"], crs=4326)
    )
    df = df.drop(columns=["lon", "lat"])

    # write to disk
    os.makedirs(os.path.dirname(ibtracs_parquet_path), exist_ok=True)
    df.to_parquet(ibtracs_parquet_path)
