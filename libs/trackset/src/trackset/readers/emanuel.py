"""
Raw ingest for Kerry Emanuel / WindRiskTech synthetic tracks.

Distributed on request as MATLAB files of ragged, zero-padded arrays, one
file per input basin ("AP", "IO", "SH", "WP" -- re-tagged here with the
STORM basin codes). Unlike CHAZ, radius to maximum winds and minimum
pressure are provided.

Ported from https://github.com/thomas-fred/emanuel-track-parser. The
end-to-end recipe (see that repo's workflow) is: ingest each basin file with
:func:`read_emanuel_mat`, concatenate, then calibrate per-basin annual
frequencies with :func:`trackset.frequency.calibrate_frequency` (against a
baseline-epoch file set) and build a TrackSet with
:func:`trackset.frequency.finalise`.

Requires scipy (``pip install nismod-trackset[emanuel]``) to read from file.
"""

import logging
from os import PathLike
from typing import Mapping, Union

import numpy as np
import pandas as pd

from ..basins import tag_basin
from ..ops import saffir_simpson_category
from ..units import MS_PER_KNOT

logger = logging.getLogger(__name__)

#: MATLAB variables holding per-(track, timestep) data
_TRACK_ARRAYS = (
    "daystore",
    "monthstore",
    "hourstore",
    "vstore",
    "rmstore",
    "pstore",
    "longstore",
    "latstore",
)


def parse_emanuel_arrays(x: Mapping[str, np.ndarray], basins=None) -> pd.DataFrame:
    """
    Convert loaded MATLAB track arrays to a dense table, tagged with basin.

    ``x`` maps MATLAB variable names to arrays of shape (n_tracks,
    max_n_timesteps) -- ``vstore`` winds in knots, ``rmstore`` radii in km,
    ``pstore`` pressures in hPa, plus per-track ``yearstore``. Padding is
    identified by ``daystore == 0``; each track's first ``n`` valid
    positions are read, where ``n`` is its count of non-padding entries
    (matching the upstream parser).
    """

    n_tracks, _ = x["daystore"].shape

    per_track = []
    for track_idx in range(n_tracks):
        n_obs = int((x["daystore"][track_idx, :] > 0).sum())
        if n_obs == 0:
            continue
        timesteps = np.arange(n_obs)
        year = int(x["yearstore"][track_idx])
        values = {
            name: x[name][track_idx, :n_obs].astype(float) for name in _TRACK_ARRAYS
        }

        wind_ms = values["vstore"] * MS_PER_KNOT
        per_track.append(
            pd.DataFrame(
                {
                    "time_utc": pd.to_datetime(
                        {
                            "year": np.full(n_obs, year),
                            "month": values["monthstore"].astype(int),
                            "day": values["daystore"].astype(int),
                            "hour": values["hourstore"].astype(int),
                        }
                    ),
                    "source_year": year,
                    "tc_number": track_idx,
                    "timestep": timesteps,
                    "sample": 0,
                    "ss_category": saffir_simpson_category(wind_ms),
                    "max_wind_speed_ms": wind_ms,
                    "radius_to_max_winds_km": values["rmstore"],
                    "min_pressure_hpa": values["pstore"],
                    "longitude_deg": values["longstore"],
                    "latitude_deg": values["latstore"],
                }
            )
        )

    df = pd.concat(per_track, ignore_index=True)

    df = tag_basin(df, lon=df["longitude_deg"], lat=df["latitude_deg"], basins=basins)

    df["track_id"] = (
        df["basin_id"]
        + df["source_year"].map(lambda x: f"_{x:04d}")
        + df["tc_number"].map(lambda x: f"_{x:04d}")
    )

    df = df.set_index("time_utc")
    df = df.sort_values(["source_year", "tc_number", "timestep"])
    df = (
        df.reset_index()
        .drop_duplicates(subset=["time_utc", "track_id"])
        .set_index("time_utc")
    )

    logger.info("Parsed %d observations from %d tracks", len(df), n_tracks)

    return df


def read_emanuel_mat(path: Union[str, PathLike], basins=None) -> pd.DataFrame:
    """
    Read one basin's MATLAB track file (see :func:`parse_emanuel_arrays`).
    """

    from scipy.io import loadmat

    return parse_emanuel_arrays(loadmat(path, squeeze_me=True), basins=basins)
