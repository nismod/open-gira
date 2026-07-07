import geopandas as gpd
import pandas as pd
import pytest


@pytest.fixture
def valid_track_points() -> gpd.GeoDataFrame:
    """
    Two tracks: 'a' heads east along the equator (4 points), 'b' sits far
    away to the north-east (2 points).
    """

    lon = [0.0, 1.0, 2.0, 3.0, 50.0, 51.0]
    lat = [0.0, 0.0, 0.0, 0.0, 30.0, 30.0]
    df = pd.DataFrame(
        {
            "track_id": ["a"] * 4 + ["b"] * 2,
            "timestep": [0, 1, 2, 3, 0, 1],
            "max_wind_speed_ms": [20.0, 35.0, 45.0, 30.0, 60.0, 72.0],
            "min_pressure_hpa": [990.0, 980.0, 970.0, 985.0, 940.0, 920.0],
            "radius_to_max_winds_km": [30.0, 25.0, 20.0, 28.0, 15.0, 12.0],
        },
        index=pd.date_range("2000-01-01", periods=6, freq="3h"),
    )
    return gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(lon, lat, crs=4326))
