"""
Operations on track point GeoDataFrames.
"""

from typing import Iterable

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry.base import BaseGeometry


#: Saffir-Simpson category lower bounds, 1-minute sustained wind speed in m/s.
#: Below 18 m/s is a tropical disturbance (-1), 18-33 a tropical storm (0).
SAFFIR_SIMPSON_LOWER_BOUNDS_MS = (18.0, 33.0, 43.0, 50.0, 58.0, 70.0)


def saffir_simpson_category(wind_speed_ms) -> np.ndarray:
    """
    Classify 1-minute sustained wind speeds (m/s) on the Saffir-Simpson scale.

    Vectorised. Returns a float array: -1 for tropical disturbance, 0 for
    tropical storm, 1-5 for hurricane categories, NaN for NaN input.

    Raises ValueError for negative wind speeds.
    """

    wind_speed_ms = np.asarray(wind_speed_ms, dtype=float)
    if (wind_speed_ms < 0).any():
        raise ValueError("wind speeds must be positive-valued")

    category = np.digitize(wind_speed_ms, SAFFIR_SIMPSON_LOWER_BOUNDS_MS) - 1.0
    return np.where(np.isnan(wind_speed_ms), np.nan, category)


def wrap_longitude(lon):
    """
    Wrap longitudes (degrees, any range, e.g. 0-360) into [-180, 180).
    """

    return ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0


def synthetic_time_index(
    df: pd.DataFrame, freq: str, start: str = "2000-01-01"
) -> pd.DataFrame:
    """
    Return a copy of ``df`` indexed by an artificial ``DatetimeIndex``.

    Synthetic track sets report a timestep number but no calendar time. Wind
    field estimation interpolates tracks and measures translation speed, which
    is easiest against a temporal index, so we fabricate one: each track's
    points are stamped from ``start`` at the set's native reporting ``freq``.

    ``df`` must have a ``track_id`` column with each track's rows stored
    contiguously and in timestep order.
    """

    lengths = df.groupby("track_id", sort=False).size()
    per_track = [
        pd.date_range(start=start, periods=length, freq=freq).values
        for length in lengths
    ]
    out = df.copy()
    out.index = pd.DatetimeIndex(np.concatenate(per_track))
    return out


def drop_duplicate_points(df: pd.DataFrame) -> pd.DataFrame:
    """
    Drop duplicated (track_id, timestep) pairs, keeping the first occurrence.

    Some synthetic sets (e.g. STORM) contain duplicate track points.
    """

    return df[~df.duplicated(subset=["track_id", "timestep"], keep="first")]


def subset_by_geometry(
    tracks: gpd.GeoDataFrame,
    geometry: BaseGeometry,
    buffer_deg: float = 0.0,
) -> gpd.GeoDataFrame:
    """
    Subset tracks to those passing within ``buffer_deg`` of ``geometry``,
    keeping each track's points from first arrival to last departure
    (including any excursions outside the area in between).

    Tracks with fewer than two points in the area of interest are dropped, as
    are tracks whose kept points are not contiguous in timestep.
    """

    area_of_interest = geometry.buffer(buffer_deg) if buffer_deg else geometry
    in_area = tracks[tracks.intersects(area_of_interest)]

    kept: list[gpd.GeoDataFrame] = []
    for track_id, points_in_area in in_area.groupby("track_id", sort=False):
        if len(points_in_area) < 2:
            continue
        track = tracks[tracks["track_id"] == track_id]
        arrival, departure = (
            points_in_area["timestep"].min(),
            points_in_area["timestep"].max(),
        )
        passage = track[track["timestep"].between(arrival, departure)]
        kept.append(passage)

    if not kept:
        return tracks.iloc[0:0]

    return pd.concat(kept)


def concat_tracks(frames: Iterable[gpd.GeoDataFrame]) -> gpd.GeoDataFrame:
    """
    Concatenate track point frames, requiring track_ids to be disjoint.
    """

    frames = list(frames)
    combined = pd.concat(frames)
    n_tracks = sum(frame["track_id"].nunique() for frame in frames)
    if combined["track_id"].nunique() != n_tracks:
        raise ValueError("track_id values collide between frames")
    return combined
