import numpy as np
import pandas as pd
import pytest
import shapely

from trackset import ops


class TestSaffirSimpsonCategory:
    def test_boundaries(self):
        # lower bound of each class: disturbance, storm, then categories 1-5
        wind_speeds = [0.0, 17.9, 18.0, 33.0, 43.0, 50.0, 58.0, 70.0, 100.0]
        expected = [-1, -1, 0, 1, 2, 3, 4, 5, 5]
        assert ops.saffir_simpson_category(wind_speeds).tolist() == expected

    def test_nan_passes_through(self):
        result = ops.saffir_simpson_category([np.nan, 20.0])
        assert np.isnan(result[0])
        assert result[1] == 0

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            ops.saffir_simpson_category([-1.0])


class TestWrapLongitude:
    def test_wrapping(self):
        assert ops.wrap_longitude(200.0) == -160.0
        assert ops.wrap_longitude(170.0) == 170.0
        assert ops.wrap_longitude(-180.0) == -180.0
        assert ops.wrap_longitude(180.0) == -180.0
        assert ops.wrap_longitude(360.0) == 0.0

    def test_vectorised(self):
        result = ops.wrap_longitude([350.0, 10.0])
        np.testing.assert_allclose(result, [-10.0, 10.0])


class TestSyntheticTimeIndex:
    def test_per_track_index(self):
        df = pd.DataFrame(
            {
                "track_id": ["a", "a", "a", "b", "b"],
                "timestep": [0, 1, 2, 0, 1],
            }
        )
        out = ops.synthetic_time_index(df, freq="3h")

        # each track's clock restarts and ticks at the native frequency
        assert out.index[0] == pd.Timestamp("2000-01-01 00:00")
        assert out.index[2] == pd.Timestamp("2000-01-01 06:00")
        assert out.index[3] == pd.Timestamp("2000-01-01 00:00")
        assert out.index[4] == pd.Timestamp("2000-01-01 03:00")


def test_drop_duplicate_points():
    df = pd.DataFrame(
        {
            "track_id": ["a", "a", "a"],
            "timestep": [0, 1, 1],
            "value": [1, 2, 3],
        }
    )
    out = ops.drop_duplicate_points(df)
    assert out["value"].tolist() == [1, 2]


class TestSubsetByGeometry:
    def test_track_through_area_kept_in_full_passage(self, valid_track_points):
        # box covers track 'a' points at timesteps 1 and 3, but not 2:
        # the passage from first arrival to last departure keeps 2 as well
        area = shapely.union(
            shapely.box(0.5, -0.5, 1.5, 0.5), shapely.box(2.5, -0.5, 3.5, 0.5)
        )
        out = ops.subset_by_geometry(valid_track_points, area)
        assert out["track_id"].unique().tolist() == ["a"]
        assert out["timestep"].tolist() == [1, 2, 3]

    def test_single_point_track_dropped(self, valid_track_points):
        # only one point of 'a' inside: dropped
        area = shapely.box(0.5, -0.5, 1.5, 0.5)
        out = ops.subset_by_geometry(valid_track_points, area)
        assert len(out) == 0

    def test_buffer_widens_selection(self, valid_track_points):
        area = shapely.box(0.6, 0.1, 2.4, 0.5)  # near but not touching 'a'
        assert len(ops.subset_by_geometry(valid_track_points, area)) == 0
        out = ops.subset_by_geometry(valid_track_points, area, buffer_deg=0.2)
        assert set(out["track_id"]) == {"a"}

    def test_no_matches_returns_empty_with_columns(self, valid_track_points):
        out = ops.subset_by_geometry(
            valid_track_points, shapely.box(-90, -50, -80, -40)
        )
        assert len(out) == 0
        assert list(out.columns) == list(valid_track_points.columns)


def test_concat_tracks_rejects_track_id_collision(valid_track_points):
    with pytest.raises(ValueError):
        ops.concat_tracks([valid_track_points, valid_track_points])
