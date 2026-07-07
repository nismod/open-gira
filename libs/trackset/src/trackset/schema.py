"""
The canonical track point schema.

A track set is a GeoDataFrame of track points with:

- a ``DatetimeIndex`` (real timestamps for observed data, synthetic ones for
  synthetic track sets -- see ``TrackSet.synthetic_time``)
- point geometry in EPSG:4326, longitudes in [-180, 180)
- one row per (track, timestep)

Required and optional columns are listed in ``REQUIRED`` and ``OPTIONAL``,
mapping column name to a dtype family ("string", "integer", "floating" or
"bool", checked with the ``pandas.api.types`` predicates).
"""

import geopandas as gpd
import numpy as np
import pandas as pd


#: Columns every track set must provide, mapped to expected dtype family.
REQUIRED: dict[str, str] = {
    "track_id": "string",  # unique per track within a set, e.g. "NA_0_1979_3"
    "timestep": "integer",  # contiguous and increasing within a track
    "max_wind_speed_ms": "floating",  # 1-minute sustained, at 10m, in m s-1
    "min_pressure_hpa": "floating",  # eye pressure in hPa / mbar
    "radius_to_max_winds_km": "floating",  # eye to max wind speed distance
}

#: Columns a track set may provide (readers preserve them when the source has
#: them; consumers must not require them).
OPTIONAL: dict[str, str] = {
    "year": "integer",
    "month": "integer",
    "tc_number": "integer",  # cardinal number of storm within year
    "basin_id": "string",  # "EP" | "NA" | "NI" | "SI" | "SP" | "WP"
    "sample": "integer",  # sample / ensemble member within a track set
    "category": "integer",  # Saffir-Simpson, -1 (disturbance) to 5
    "landfall": "bool",
    "distance_to_land_km": "floating",
    "name": "string",  # observed storms may be named
}

_DTYPE_CHECKS = {
    "string": pd.api.types.is_string_dtype,
    "integer": pd.api.types.is_integer_dtype,
    "floating": pd.api.types.is_float_dtype,
    "bool": pd.api.types.is_bool_dtype,
}


class SchemaError(ValueError):
    """A track set GeoDataFrame does not conform to the canonical schema."""


def problems(gdf: gpd.GeoDataFrame) -> list[str]:
    """
    Check a GeoDataFrame against the canonical track point schema.

    Returns a list of human-readable problems; empty if the frame conforms.
    """

    found: list[str] = []

    if not isinstance(gdf, gpd.GeoDataFrame):
        return [f"expected GeoDataFrame, got {type(gdf)}"]

    if gdf.crs is None or gdf.crs.to_epsg() != 4326:
        found.append(f"CRS must be EPSG:4326, got {gdf.crs}")

    if not isinstance(gdf.index, pd.DatetimeIndex):
        found.append(f"index must be DatetimeIndex, got {type(gdf.index).__name__}")

    for name, family in REQUIRED.items():
        if name not in gdf.columns:
            found.append(f"missing required column: {name}")
        elif not _DTYPE_CHECKS[family](gdf[name].dtype):
            found.append(f"column {name} should be {family}, got {gdf[name].dtype}")
        elif gdf[name].isna().any():
            found.append(f"required column {name} contains null values")

    for name, family in OPTIONAL.items():
        if name in gdf.columns and not _DTYPE_CHECKS[family](gdf[name].dtype):
            found.append(f"column {name} should be {family}, got {gdf[name].dtype}")

    # geometry and per-track checks only make sense on non-empty, structurally
    # sound frames
    if len(gdf) == 0 or found:
        return found

    if not (gdf.geometry.geom_type == "Point").all():
        found.append("geometry must be points")
    else:
        lon = gdf.geometry.x
        if (lon < -180).any() or (lon >= 180).any():
            found.append("longitudes must lie in [-180, 180)")

    for track_id, track in gdf.groupby("track_id", sort=False):
        steps = track["timestep"].to_numpy()
        if not (np.diff(steps) == 1).all():
            found.append(
                f"track {track_id}: timestep must be contiguous and increasing"
            )

    return found


def validate(gdf: gpd.GeoDataFrame) -> None:
    """
    Raise ``SchemaError`` listing all problems if ``gdf`` does not conform to
    the canonical track point schema.
    """

    found = problems(gdf)
    if found:
        raise SchemaError(
            "track set does not conform to schema:\n- " + "\n- ".join(found)
        )
