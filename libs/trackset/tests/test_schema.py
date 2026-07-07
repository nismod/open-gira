import pandas as pd
import pytest

from trackset import SchemaError, problems, validate


def test_valid_frame_passes(valid_track_points):
    validate(valid_track_points)


def test_empty_frame_with_columns_passes(valid_track_points):
    validate(valid_track_points.iloc[0:0])


def test_not_a_geodataframe():
    assert problems(pd.DataFrame({"track_id": []})) != []


def test_missing_required_column(valid_track_points):
    broken = valid_track_points.drop(columns=["radius_to_max_winds_km"])
    assert any("radius_to_max_winds_km" in p for p in problems(broken))


def test_null_in_required_column(valid_track_points):
    broken = valid_track_points.copy()
    broken.loc[broken.index[0], "max_wind_speed_ms"] = None
    assert any("null" in p for p in problems(broken))


def test_wrong_crs(valid_track_points):
    broken = valid_track_points.to_crs(3857)
    assert any("EPSG:4326" in p for p in problems(broken))


def test_longitude_out_of_range(valid_track_points):
    broken = valid_track_points.copy()
    broken.geometry = broken.geometry.translate(xoff=340.0)
    assert any("longitude" in p for p in problems(broken))


def test_non_datetime_index(valid_track_points):
    broken = valid_track_points.reset_index(drop=True)
    assert any("DatetimeIndex" in p for p in problems(broken))


def test_non_contiguous_timestep(valid_track_points):
    broken = valid_track_points.copy()
    broken.loc[broken.index[2], "timestep"] = 5
    assert any("contiguous" in p for p in problems(broken))


def test_validate_raises_with_all_problems(valid_track_points):
    broken = valid_track_points.drop(columns=["track_id", "timestep"])
    with pytest.raises(SchemaError) as excinfo:
        validate(broken)
    assert "track_id" in str(excinfo.value)
    assert "timestep" in str(excinfo.value)
