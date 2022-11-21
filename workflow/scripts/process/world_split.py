"""Splits the world into boxes of equal length and height
"""
import json
import sys

import geopandas
import numpy
from shapely.geometry import Point


def create_global_grid(box_width_height):
    assert 180 % box_width_height == 0  # ensure divisibility

    lat_centroids = numpy.arange(90 - box_width_height / 2, -90, -box_width_height)
    lon_centroids = numpy.arange(-180 + box_width_height / 2, 180, box_width_height)

    points = []
    for lat in lat_centroids:  # left to right and down
        for lon in lon_centroids:
            points.append(Point(lon, lat))
    points_gdf = geopandas.GeoDataFrame({"geometry": points})
    grid_geoms = points_gdf.buffer(box_width_height / 2, cap_style=3)

    grid = geopandas.GeoDataFrame(
        {
            "geometry": grid_geoms,
            "box_id": [f"box_{i}" for i in range(len(grid_geoms))],
        },
        crs="EPSG:4326",
    )
    ncols = len(lon_centroids)
    nrows = len(lat_centroids)
    return grid, ncols, nrows


if __name__ == "__main__":
    try:
        admin_data_path = snakemake.input["admin_data"]  # type: ignore
        box_width_height = snakemake.config["box_deg"]  # type: ignore
        output_dir = snakemake.config["output_dir"]  # type: ignore
        global_metadata_path = snakemake.output["global_metadata"]  # type: ignore
        global_boxes_path = snakemake.output["global_boxes"]  # type: ignore
    except:
        admin_data_path = sys.argv[1]
        output_dir = sys.argv[2]
        box_width_height = sys.argv[3]
        global_metadata_path = sys.argv[4]
        global_boxes_path = sys.argv[5]

    box_width_height = float(box_width_height)
    grid, ncols, nrows = create_global_grid(box_width_height)

    countries = geopandas.read_parquet(admin_data_path) \
        .drop(["NAME_0"], axis="columns") \
        .rename({"GID_0": "code"}, axis="columns")

    # e.g. "box_284" -> ['GBR', 'FRA']
    countries_by_box: dict[str, list[str]] = {}
    for box_id, box_countries in geopandas.sjoin(countries, grid).groupby('box_id'):
        countries_by_box[box_id] = list(box_countries.code)

    # post-processing to match Max's previous output
    # 1) assign None to empty boxes
    for box in grid.itertuples():
        if box.box_id not in countries_by_box:
            countries_by_box[box.box_id] = None
    # 2) sort the elements by their box index
    countries_by_box = dict(sorted(countries_by_box.items(), key=lambda item: int(item[0].split("_")[-1])))

    with open(global_metadata_path, "w") as filejson:
        lon_min, lat_min, lon_max, lat_max = grid.bounds.values[0]
        info = {
            "boxlen": box_width_height,
            "lon_min": lon_min,
            "lat_min": lat_min,
            "lon_max": lon_max,
            "lat_max": lat_max,
            "num_cols": ncols,
            "num_rows": nrows,
            "tot_boxes": int(ncols * nrows),
            "box_country_dict": countries_by_box,
        }
        json.dump(info, filejson, indent=2, sort_keys=True)

    grid.to_parquet(global_boxes_path)
