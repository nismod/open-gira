import pandas as pd

from trackset.basins import basin_polygons, tag_basin


def test_basin_polygons():
    basins = basin_polygons()
    assert set(basins["basin_id"]) == {"EP", "NA", "NI", "SI", "SP", "WP"}
    assert basins.crs.to_epsg() == 4326


def test_tag_basin_signed_and_positive_longitudes():
    df = pd.DataFrame({"point": ["caribbean", "fiji", "fiji_positive", "tokyo"]})
    # Fiji given in both signed (-178) and 0-360 (182) conventions
    lon = [-60.0, -178.0, 182.0, 140.0]
    lat = [20.0, -18.0, -18.0, 30.0]

    tagged = tag_basin(df, lon=lon, lat=lat)

    assert tagged.set_index("point")["basin_id"].to_dict() == {
        "caribbean": "NA",
        "fiji": "SP",
        "fiji_positive": "SP",
        "tokyo": "WP",
    }
    assert "geometry" not in tagged.columns


def test_tag_basin_drops_points_outside_basins():
    df = pd.DataFrame({"point": ["north_sea", "gulf_of_mexico"]})
    tagged = tag_basin(df, lon=[3.0, -90.0], lat=[56.0, 25.0])
    # 3E lies in the 0-10E gap between the NA and NI/SI basins: dropped
    assert tagged["point"].tolist() == ["gulf_of_mexico"]
    assert tagged["basin_id"].tolist() == ["NA"]
