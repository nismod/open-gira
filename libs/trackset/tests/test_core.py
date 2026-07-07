import pytest
import shapely

from trackset import SchemaError, TrackSet


def make_trackset(gdf, **overrides) -> TrackSet:
    kwargs = dict(
        data=gdf,
        source="TEST",
        years=10.0,
        synthetic_time=True,
        attributes={"gcm": "TEST-GCM", "epoch": 2050},
    )
    kwargs.update(overrides)
    return TrackSet(**kwargs)


def test_construction_validates(valid_track_points):
    broken = valid_track_points.drop(columns=["timestep"])
    with pytest.raises(SchemaError):
        make_trackset(broken)


def test_years_must_be_positive(valid_track_points):
    with pytest.raises(ValueError):
        make_trackset(valid_track_points, years=0.0)


def test_annual_frequency(valid_track_points):
    ts = make_trackset(valid_track_points)  # 2 tracks over 10 years
    assert ts.n_tracks == 2
    assert ts.annual_frequency == pytest.approx(0.2)


def test_subset_preserves_metadata(valid_track_points):
    ts = make_trackset(valid_track_points)
    subset = ts.subset(shapely.box(-1, -1, 4, 1))  # covers track 'a' only
    assert subset.n_tracks == 1
    assert subset.years == ts.years
    assert subset.source == ts.source


def test_parquet_round_trip(tmp_path, valid_track_points):
    ts = make_trackset(valid_track_points)
    path = tmp_path / "tracks.geoparquet"
    ts.to_parquet(path)

    recovered = TrackSet.read_parquet(path)

    assert recovered.source == ts.source
    assert recovered.years == ts.years
    assert recovered.synthetic_time == ts.synthetic_time
    assert recovered.wind_averaging_period == ts.wind_averaging_period
    assert recovered.attributes == ts.attributes
    assert recovered.data.equals(ts.data)


def test_read_parquet_rejects_plain_geoparquet(tmp_path, valid_track_points):
    path = tmp_path / "plain.geoparquet"
    valid_track_points.to_parquet(path)
    with pytest.raises(ValueError, match="metadata"):
        TrackSet.read_parquet(path)
