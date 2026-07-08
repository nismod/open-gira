import geopandas as gpd
import numpy as np
import pandas as pd
import pytest

from trackset import TrackSet, frequency


def make_synthetic_tracks(
    n_tracks_per_basin: dict[str, int], points_per_track: int = 3
) -> pd.DataFrame:
    """
    Fabricate a raw synthetic track table (as produced by the CHAZ/Emanuel
    ingest readers) with the given number of tracks per basin.
    """

    frames = []
    basin_centre_lon = {"NA": -60.0, "WP": 140.0, "SP": 175.0}
    for basin_id, n_tracks in n_tracks_per_basin.items():
        for i in range(n_tracks):
            source_year = 2000 + (i % 20)
            frames.append(
                pd.DataFrame(
                    {
                        "track_id": f"{basin_id}_{source_year:04d}_{i:04d}",
                        "source_year": source_year,
                        "timestep": range(points_per_track),
                        "ss_category": 1.0,
                        "max_wind_speed_ms": 35.0,
                        "min_pressure_hpa": 980.0,
                        "radius_to_max_winds_km": 30.0,
                        "basin_id": basin_id,
                        "longitude_deg": basin_centre_lon[basin_id] + 0.1 * i,
                        "latitude_deg": 15.0 if basin_id != "SP" else -15.0,
                    },
                    index=pd.date_range(
                        f"{source_year}-06-01", periods=points_per_track, freq="3h"
                    ),
                )
            )
    return pd.concat(frames)


class TestTracksPerYear:
    def test_counts_unique_tracks_over_span(self):
        # 40 NA tracks and 20 WP tracks with source years spanning 2000-2019
        tracks = make_synthetic_tracks({"NA": 40, "WP": 20})
        result = frequency.tracks_per_year(tracks)
        assert result["NA"] == pytest.approx(2.0)
        assert result["WP"] == pytest.approx(1.0)


def test_relative_frequency():
    baseline = make_synthetic_tracks({"NA": 20, "WP": 20})
    target = make_synthetic_tracks({"NA": 40, "WP": 20})
    change = frequency.relative_frequency(baseline, target)
    assert change["NA"] == pytest.approx(2.0)
    assert change["WP"] == pytest.approx(1.0)


def test_observed_frequency():
    # observed-style tracks: geometry + year columns, agency basin labels to
    # be replaced by polygon tagging
    n_years = (
        frequency.DEFAULT_OBSERVATION_WINDOW[1]
        - frequency.DEFAULT_OBSERVATION_WINDOW[0]
        + 1
    )
    points = []
    for i in range(n_years):
        year = frequency.DEFAULT_OBSERVATION_WINDOW[0] + i
        points.append({"track_id": f"storm_{year}", "year": year, "lon": -60.0})
    # one storm outside the window: ignored
    points.append({"track_id": "storm_old", "year": 1980, "lon": -60.0})
    df = pd.DataFrame(points)
    tracks = gpd.GeoDataFrame(
        df[["track_id", "year"]],
        geometry=gpd.points_from_xy(df["lon"], [20.0] * len(df)),
        crs=4326,
    )

    observed = frequency.observed_frequency(tracks)

    # one storm per year in the North Atlantic
    assert observed["NA"] == pytest.approx(1.0)


class TestCalibrateFrequency:
    def test_calibration(self):
        tracks = make_synthetic_tracks({"NA": 100, "WP": 50})
        target = pd.Series({"NA": 10.0, "WP": 2.0}, name="tc_per_year")

        calibrated, duration = frequency.calibrate_frequency(
            tracks, target, rng=np.random.default_rng(0)
        )

        # NA can represent 100/10 = 10 years, WP 50/2 = 25: truncate to 10
        assert duration == 10.0
        assert (calibrated["year"] < duration).all()
        # sample is the millennium chunk of the assigned year
        assert (calibrated["sample"] == calibrated["year"] // 1000).all()
        # tc_number unique per (year, track)
        per_year = calibrated.loc[
            :, ["year", "tc_number", "track_id"]
        ].drop_duplicates()
        assert not per_year.duplicated(subset=["year", "tc_number"]).any()
        # all original columns preserved
        assert "max_wind_speed_ms" in calibrated.columns

    def test_deterministic_with_seed(self):
        tracks = make_synthetic_tracks({"NA": 30, "WP": 30})
        target = pd.Series({"NA": 3.0, "WP": 3.0}, name="tc_per_year")

        first, _ = frequency.calibrate_frequency(
            tracks, target, rng=np.random.default_rng(42)
        )
        second, _ = frequency.calibrate_frequency(
            tracks, target, rng=np.random.default_rng(42)
        )
        pd.testing.assert_frame_equal(first, second)

    def test_rejects_basin_without_tracks(self):
        tracks = make_synthetic_tracks({"NA": 10})
        target = pd.Series({"NA": 1.0, "WP": 1.0}, name="tc_per_year")
        with pytest.raises(ValueError, match="no tracks"):
            frequency.calibrate_frequency(tracks, target, rng=np.random.default_rng(0))


class TestFinalise:
    def calibrated(self) -> tuple[pd.DataFrame, float]:
        tracks = make_synthetic_tracks({"NA": 50, "WP": 50})
        target = pd.Series({"NA": 5.0, "WP": 5.0}, name="tc_per_year")
        return frequency.calibrate_frequency(
            tracks, target, rng=np.random.default_rng(1)
        )

    def test_finalise_builds_valid_trackset(self):
        calibrated, duration = self.calibrated()
        ts = frequency.finalise(
            calibrated, source="TEST_scenario", years=duration, attributes={"gcm": "X"}
        )
        assert isinstance(ts, TrackSet)
        assert ts.years == duration
        assert ts.synthetic_time
        assert ts.attributes == {"gcm": "X"}
        # ss_category renamed and integral
        assert ts.data["category"].dtype.kind == "i"

    def test_finalise_drops_points_and_broken_tracks(self):
        calibrated, duration = self.calibrated()

        # null out a mid-track pressure: the point is invalid, and its
        # removal leaves the track with a timestep gap, so the track goes
        a_track = calibrated["track_id"].iloc[0]
        mid_point = (calibrated["track_id"] == a_track) & (calibrated["timestep"] == 1)
        calibrated.loc[mid_point, "min_pressure_hpa"] = np.nan

        ts = frequency.finalise(calibrated, source="TEST", years=duration)
        assert a_track not in set(ts.data["track_id"])

    def test_finalise_requires_schema_columns(self):
        calibrated, duration = self.calibrated()
        with pytest.raises(ValueError, match="required"):
            frequency.finalise(
                calibrated.drop(columns=["min_pressure_hpa"]),
                source="TEST",
                years=duration,
            )
