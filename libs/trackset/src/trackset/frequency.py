"""
Per-basin annual frequency calculation and calibration of synthetic track
sets against observations.

Synthetic track generators do not necessarily produce storms at realistic
rates. The approach here (ported from thomas-fred/chaz-track-parser and
thomas-fred/emanuel-track-parser, where it was developed to calibrate CHAZ
and Emanuel track sets):

1. measure observed tracks-per-year per basin from IBTrACS
   (:func:`observed_frequency`);
2. for future epochs, scale by the *relative* change in frequency between
   the synthetic set's own baseline and target epochs
   (:func:`relative_frequency`) -- trusting the generator for the climate
   change signal but not the absolute rate;
3. re-assign each synthetic track to a random year drawn from the longest
   duration consistent with the target frequency, giving a set with a known,
   physically-anchored time span (:func:`calibrate_frequency`).

Unlike the upstream implementations, :func:`calibrate_frequency` requires an
explicitly seeded random generator: the year re-assignment is stochastic,
and an unseeded draw makes the output impossible to reproduce from source.
"""

import logging
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd

from . import ops
from .basins import tag_basin
from .core import TrackSet
from .schema import OPTIONAL, REQUIRED

logger = logging.getLogger(__name__)

#: Observation window for IBTrACS-derived frequencies. Data outside this
#: window shows a markedly different frequency.
DEFAULT_OBSERVATION_WINDOW = (2002, 2023)


def tracks_per_year(tracks: pd.DataFrame, year_col: str = "source_year") -> pd.Series:
    """
    Mean annual track count per basin, over the year span present in
    ``tracks`` (which must have ``track_id``, ``basin_id`` and ``year_col``
    columns).
    """

    duration = tracks[year_col].max() - tracks[year_col].min() + 1
    return (
        tracks.loc[:, ["basin_id", "track_id"]]
        .groupby(["basin_id"])
        .nunique()
        .rename(columns={"track_id": "tc_per_year"})
        .loc[:, "tc_per_year"]
        / duration
    )


def relative_frequency(
    baseline_tracks: pd.DataFrame,
    target_tracks: pd.DataFrame,
    year_col: str = "source_year",
) -> pd.Series:
    """
    Relative change in per-basin annual track frequency between two epochs of
    the same synthetic set (target / baseline).
    """

    return tracks_per_year(target_tracks, year_col) / tracks_per_year(
        baseline_tracks, year_col
    )


def observed_frequency(
    tracks: gpd.GeoDataFrame,
    window: tuple[int, int] = DEFAULT_OBSERVATION_WINDOW,
    basins: gpd.GeoDataFrame | None = None,
) -> pd.Series:
    """
    Observed annual track count per basin from a set of observed tracks
    (typically ``read_ibtracs(...).data``), within an inclusive year window.

    Points are (re-)tagged against the basin polygons, replacing any
    reporting-agency basin labels, for consistency with the tagging applied
    to synthetic tracks.
    """

    start_year, end_year = window
    tagged = tag_basin(
        tracks.drop(columns=["basin_id", "geometry"], errors="ignore"),
        lon=tracks.geometry.x,
        lat=tracks.geometry.y,
        basins=basins,
    )
    in_window = tagged[(tagged.year >= start_year) & (tagged.year <= end_year)]
    duration = end_year - start_year + 1
    return (
        in_window.loc[:, ["basin_id", "track_id"]]
        .groupby(["basin_id"])
        .nunique()
        .rename(columns={"track_id": "tc_per_year"})
        .loc[:, "tc_per_year"]
        / duration
    )


def calibrate_frequency(
    tracks: pd.DataFrame,
    target_tracks_per_year: pd.Series,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, float]:
    """
    Re-assign synthetic tracks to years such that per-basin annual
    frequencies match ``target_tracks_per_year``.

    Given the target frequency and the available track count, each basin can
    represent some duration (count / frequency); tracks are assigned uniform
    random years within their basin's duration, and the result is truncated
    to the shortest duration across basins so that every represented year
    has global coverage.

    Args:
        tracks: Synthetic track points with ``track_id``, ``basin_id`` and
            ``timestep`` columns and a DatetimeIndex. All other columns are
            preserved.
        target_tracks_per_year: Target mean annual track count, indexed by
            ``basin_id`` (e.g. from :func:`observed_frequency`, optionally
            scaled by :func:`relative_frequency`).
        rng: Seeded random generator for the year assignment, e.g.
            ``np.random.default_rng(0)``. Reruns with the same seed and
            inputs give identical output.

    Returns:
        (calibrated tracks, duration): the tracks with re-assigned ``year``,
        per-year ``tc_number``, and millennium-chunk ``sample`` columns; and
        the number of years the calibrated set represents (the ``years``
        for a :class:`TrackSet` built from it).
    """

    frequency = target_tracks_per_year.to_frame(name="tc_per_year")

    absent = set(frequency.index) - set(tracks["basin_id"].unique())
    if absent:
        raise ValueError(
            f"target frequencies given for basins with no tracks: {absent}"
        )

    # given the desired average track frequency, how many years might we represent?
    frequency["track_count"] = tracks.groupby("basin_id")["track_id"].nunique()
    frequency["duration_years"] = np.round(
        frequency["track_count"] / frequency["tc_per_year"], 0
    ).astype(int)

    track_basin = (
        tracks.loc[:, ["track_id", "basin_id"]]
        .drop_duplicates()
        .set_index("track_id", drop=True)
    )
    track_year = []
    logger.info("Resampling years by basin")
    for basin_id in frequency.index:
        basin = track_basin[track_basin.basin_id == basin_id].copy()
        basin["year"] = np.round(
            rng.random(len(basin)) * frequency.loc[basin_id, "duration_years"], 0
        ).astype(int)
        track_year.append(basin)

    # truncate to the shortest basin duration, for global coverage of every year
    duration = int(frequency["duration_years"].min())
    track_year = pd.concat(track_year).sort_values("year")
    track_year = track_year[track_year.year < duration]

    logger.info("Joining years to tracks")
    calibrated = tracks.drop(columns=["year", "tc_number"], errors="ignore").join(
        track_year.year, on="track_id", how="inner"
    )

    # label with a tc_number, unique within a given year
    track_year = track_year.reset_index().drop(columns=["basin_id"]).drop_duplicates()
    track_year["tc_number"] = track_year.groupby("year").cumcount()

    # pandas merge drops the (datetime) index; preserve it via a temporary column
    index_name = calibrated.index.name or "index"
    calibrated = (
        calibrated.rename_axis(index_name)
        .reset_index()
        .merge(track_year, on=["year", "track_id"], how="left", validate="m:1")
        .set_index(index_name)
        .rename_axis(None)
    )
    assert (calibrated["tc_number"] >= 0).all()

    calibrated = calibrated.sort_values(["year", "tc_number", "timestep"])

    # label sample as millennium chunks
    calibrated["sample"] = np.floor(calibrated["year"] / 1000).astype(int)

    # drop any duplicate track observations
    calibrated = (
        calibrated.rename_axis(index_name)
        .reset_index()
        .drop_duplicates(subset=[index_name, "track_id"])
        .set_index(index_name)
        .rename_axis(None)
    )

    return calibrated, float(duration)


def finalise(
    calibrated: pd.DataFrame,
    source: str,
    years: float,
    attributes: dict[str, Any] | None = None,
) -> TrackSet:
    """
    Build a :class:`TrackSet` from a calibrated synthetic track table.

    Expects ``longitude_deg`` / ``latitude_deg`` columns (any longitude
    range) and the calibration outputs of :func:`calibrate_frequency`.
    Points missing a required variable (e.g. an implausible inferred
    pressure) are dropped, and tracks left with non-contiguous timesteps by
    that removal are dropped entirely (the same policy open-gira applies
    when slicing tracks); both removals are logged.

    ``synthetic_time`` is set: calibration re-assigns years, so the
    source-simulation timestamps retain seasonality but not calendar meaning.
    """

    df = calibrated.rename(columns={"ss_category": "category"})

    keep = [c for c in [*REQUIRED, *OPTIONAL] if c in df.columns]
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"calibrated tracks missing required columns: {missing}")

    geometry = gpd.points_from_xy(
        ops.wrap_longitude(df["longitude_deg"]), df["latitude_deg"], crs=4326
    )
    df = gpd.GeoDataFrame(df.loc[:, keep], geometry=geometry)

    valid = df[list(REQUIRED)].notna().all(axis="columns")
    if not valid.all():
        logger.info("Dropping %d points with missing required values", (~valid).sum())
        df = df[valid]

    contiguous = df.groupby("track_id", sort=False)["timestep"].transform(
        lambda steps: (steps.diff().dropna() == 1).all()
    )
    if not contiguous.all():
        n_tracks = df.loc[~contiguous, "track_id"].nunique()
        logger.info("Dropping %d tracks with non-contiguous timesteps", n_tracks)
        df = df[contiguous]

    if "category" in df.columns:
        df["category"] = df["category"].astype(int)

    return TrackSet(
        data=df,
        source=source,
        years=years,
        synthetic_time=True,
        attributes=attributes or {},
    )
