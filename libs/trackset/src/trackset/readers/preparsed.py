"""
Readers for CHAZ (Lee et al., https://doi.org/10.1002/2017MS001186) and
Kerry Emanuel / WindRiskTech track sets.

Both are distributed on request as ragged arrays (netCDF and MATLAB
respectively) and are currently pre-processed to tabular GeoParquet by
external tools (https://github.com/thomas-fred/chaz-track-parser and
https://github.com/thomas-fred/emanuel-track-parser), which also calibrate
annual track frequencies to the historical record per-basin. Absorbing those
parsers into this package is on the roadmap; until then, these readers
normalise the pre-parsed tabular output.

Because the raw formats do not state their simulated time span, ``years``
must be supplied by the caller.
"""

from os import PathLike
from typing import Union

import geopandas as gpd

from .. import ops
from ..core import TrackSet


def _normalise_preparsed(df: gpd.GeoDataFrame, source: str, years: float) -> TrackSet:
    df = df.copy()

    df["track_id"] = (
        df["sample"].map(lambda x: f"S{x:03d}")
        + df["year"].map(lambda x: f"Y{x:04d}")
        + df["tc_number"].map(lambda x: f"N{x:03d}")
    )

    lon = ops.wrap_longitude(df.geometry.x)
    df = df.set_geometry(gpd.points_from_xy(lon, df.geometry.y, crs=4326))

    df = ops.drop_duplicate_points(df)

    return TrackSet(data=df, source=source, years=years, synthetic_time=True)


def read_chaz(path: Union[str, PathLike], source: str, years: float) -> TrackSet:
    """
    Read pre-parsed CHAZ tracks from GeoParquet.

    Args:
        path: GeoParquet output of chaz-track-parser.
        source: Scenario label, e.g. "CHAZ_SSP-585_GCM-UKESM1-0-LL_epoch-2050".
        years: Simulated years this file represents.
    """

    return _normalise_preparsed(gpd.read_parquet(path), source, years)


def read_emanuel(path: Union[str, PathLike], source: str, years: float) -> TrackSet:
    """
    Read pre-parsed Kerry Emanuel tracks from GeoParquet.

    Args:
        path: GeoParquet output of emanuel-track-parser.
        source: Scenario label, e.g. "emanuel_ssp-585_gcm-ukmo6_epoch-2050".
        years: Simulated years this file represents (typically ~200).
    """

    return _normalise_preparsed(gpd.read_parquet(path), source, years)
