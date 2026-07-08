"""
Raw ingest for CHAZ synthetic tracks (Lee et al.).

https://doi.org/10.1002/2017MS001186

CHAZ is distributed on request as ragged netCDF datacubes with dimensions
(lifelength, stormID) and, for wind speed, an additional ensemble dimension.
CHAZ reports only position and maximum wind speed: the radius to maximum
winds and minimum pressure needed by wind field models are inferred with the
parametric fits in :mod:`trackset.physics`.

Ported from https://github.com/thomas-fred/chaz-track-parser. The end-to-end
recipe (see that repo's workflow) is: ingest each sample file with
:func:`read_chaz_dataset`, window to a baseline and a target epoch with
:func:`filter_epoch`, then calibrate per-basin annual frequencies with
:func:`trackset.frequency.calibrate_frequency` using observed rates scaled
by the epochs' relative frequency, and build a TrackSet with
:func:`trackset.frequency.finalise`.

Requires xarray (``pip install nismod-trackset[chaz]``) to read from file.
"""

import logging
from math import prod
from os import PathLike
from typing import Iterator, Union

import numpy as np
import pandas as pd

from ..basins import tag_basin
from ..ops import saffir_simpson_category
from ..physics import ENV_PRESSURE, p_min_holland_1980, r_max_willoughby_2004
from ..units import MS_PER_KNOT

logger = logging.getLogger(__name__)

#: Short code for use in track_id, keyed by genesis method
GENESIS_METHOD_CODE = {
    "SD": "S",  # Saturation Deficit
    "CRH": "H",  # Column Relative Humidity
}

#: The netCDF 'time' variable is days since this date
REFERENCE_DATE = "1950-01-01"


def iter_chaz_ensembles(ds, genesis_method: str, sample: int) -> Iterator[pd.DataFrame]:
    """
    Convert a CHAZ netCDF dataset to dense tables, one per ensemble member.

    Creates timestamps from the reference date and day offsets, unravels the
    (lifelength, stormID) datacube to long format, converts winds from knots
    to m/s, builds a unique ``track_id`` and drops padding (NaN) rows.

    Args:
        ds: An ``xarray.Dataset`` of CHAZ tracks.
        genesis_method: "SD" or "CRH", used in track ids.
        sample: Sample number of this file, used in track ids.
    """

    # check ordering of dimensions
    for var in ("time", "longitude", "latitude"):
        assert ds[var].dims == ("lifelength", "stormID")
    assert ds["Mwspd"].dims == ("ensembleNum", "lifelength", "stormID")

    # check we won't overflow our fixed-width strings in track_id
    assert np.logical_and(0 <= ds.ensembleNum, ds.ensembleNum < 1e2).all()
    assert np.logical_and(0 <= ds.stormID, ds.stormID < 1e5).all()

    timestamps_2d = pd.to_datetime(
        ds.time.values.ravel(), unit="D", origin=pd.Timestamp(REFERENCE_DATE)
    ).values.reshape(ds.time.values.shape)

    length = prod(timestamps_2d.shape)
    storm = np.repeat(ds.stormID.data, ds.sizes["lifelength"])
    # N.B. We transpose the input arrays to have stormID as the first dim,
    # then lifelength
    source_year = np.repeat(
        timestamps_2d.T[:, 0].astype("datetime64[Y]").astype(int) + 1970,
        ds.sizes["lifelength"],
    )
    timesteps = np.tile(range(ds.sizes["lifelength"]), len(ds.stormID.data))
    timestamps = timestamps_2d.T.reshape(length)
    longitude = ds.longitude.data.T.reshape(length)
    latitude = ds.latitude.data.T.reshape(length)

    for i in ds.ensembleNum.data:
        df = (
            pd.DataFrame(
                {
                    "time_utc": timestamps,
                    "source_year": source_year,
                    "storm": storm,
                    "sample": np.ones(length) * int(sample),
                    "ensemble": np.ones(length) * i,
                    "timestep": timesteps,
                    "longitude_deg": longitude,
                    "latitude_deg": latitude,
                    "max_wind_speed_ms": ds.Mwspd.data[i, :, :].T.reshape(length)
                    * MS_PER_KNOT,
                }
            )
            .dropna()
            .set_index("time_utc", drop=True)
            .astype({"sample": int, "ensemble": int})
        )

        df["track_id"] = (
            f"{GENESIS_METHOD_CODE[genesis_method]}"
            + f"_{int(sample):03d}"
            + df["source_year"].map(lambda x: f"_{x:04d}")
            + df["storm"].map(lambda x: f"_{x:05d}")
            + df["ensemble"].map(lambda x: f"_{x:02d}")
        )

        df = (
            df.reset_index()
            .drop_duplicates(subset=["time_utc", "track_id"])
            .set_index("time_utc")
        )

        yield df


def process_chaz_ensemble(df: pd.DataFrame, basins=None) -> pd.DataFrame:
    """
    Enrich a dense CHAZ ensemble table: tag Saffir-Simpson category and
    basin, and infer radius to maximum winds (Willoughby 2004) and minimum
    pressure (Holland 1980 profile, Vickery & Wadhera 2008 shape parameter).

    Inferred pressures outside the plausible range are NaN; see
    :func:`trackset.frequency.finalise` for the removal policy.
    """

    df = df.copy()
    df["ss_category"] = saffir_simpson_category(df["max_wind_speed_ms"])
    df = tag_basin(df, lon=df["longitude_deg"], lat=df["latitude_deg"], basins=basins)
    df["radius_to_max_winds_km"] = r_max_willoughby_2004(
        df["max_wind_speed_ms"], df["latitude_deg"]
    )
    df["min_pressure_hpa"] = p_min_holland_1980(
        df["basin_id"].map(ENV_PRESSURE),
        df["max_wind_speed_ms"],
        df["radius_to_max_winds_km"] * 1_000,
        df["latitude_deg"],
    )
    return df


def read_chaz_dataset(
    ds, genesis_method: str, sample: int, basins=None
) -> pd.DataFrame:
    """
    Read a CHAZ ``xarray.Dataset``: all ensemble members as one dense,
    enriched table (see :func:`iter_chaz_ensembles` and
    :func:`process_chaz_ensemble`). For very large files, iterate the
    generator and process/write per ensemble instead.
    """

    ensembles = []
    for chunk in iter_chaz_ensembles(ds, genesis_method, sample):
        ensembles.append(process_chaz_ensemble(chunk, basins=basins))
        logger.info("Processed ensemble %d", ensembles[-1]["ensemble"].iloc[0])
    return pd.concat(ensembles)


def read_chaz_netcdf(
    path: Union[str, PathLike], genesis_method: str, sample: int, basins=None
) -> pd.DataFrame:
    """
    Read a CHAZ netCDF file (see :func:`read_chaz_dataset`).
    """

    import xarray as xr

    return read_chaz_dataset(
        xr.open_dataset(path), genesis_method, sample, basins=basins
    )


def filter_epoch(df: pd.DataFrame, epoch: int, half_width_years: int) -> pd.DataFrame:
    """
    Subset track points to a window of ``source_year`` values (exclusive)
    around an epoch, e.g. 2050 +/- 15 years.
    """

    lower = epoch - half_width_years
    upper = epoch + half_width_years
    if lower < df.source_year.min() or df.source_year.max() < upper:
        raise ValueError(f"{epoch:d} +/- {half_width_years} is out of range")
    return df[(lower < df.source_year) & (df.source_year < upper)]
