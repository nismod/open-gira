"""
Given some geospatial data, identify which administrative area each row of the
input data occupies. For a specified set of fields, determine the row-wise
aggregate (typically sum) for every field, for every admin area.

Available administrative layers in GADM:
    - National (level 0)
    - State/province/equivalent (level 1)
    - County/district/equivalent (level 2)
    - Smaller Level 3 or 4.

Seldom does a country have boundaries for every level.

N.B. If geometries intersect with multiple admin areas, they will contribute
the aggregate of each selected field in full, to every intersected admin area.
This means there is some 'double counting' especially when aggregating to
smaller administrative subdivisions such as levels 2, 3 and 4.
"""


import json
import logging
import warnings
import re
import sys

import geopandas as gpd
import pyproj
import shapely


WGS84_EPSG = 4326


if __name__ == "__main__":

    try:
        admin_areas_path: str = snakemake.input["admin_areas"]
        slice_bounds_path: str = snakemake.input["slice_bounds"]
        data_to_aggregate_path: str = snakemake.input["data_to_aggregate"]
        column_regex: str = snakemake.params["columns_to_aggregate_regex"]
        admin_level_slug: int = snakemake.wildcards.ADMIN_SLUG
        agg_func_slug: str = snakemake.wildcards.AGG_FUNC_SLUG
        aggregated_path: str = snakemake.output["aggregated_data"]
    except NameError:
        raise ValueError("Must be run via snakemake.")

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    # N.B. geoparquet spec is not currently stable, and therefore not a suitable format for archiving
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    logging.info("Reading data to aggregate")
    data_to_aggregate: gpd.GeoDataFrame = gpd.read_parquet(data_to_aggregate_path)

    with open(slice_bounds_path, "r") as fp:
        extracts, = json.load(fp)["extracts"]
        minx, miny, maxx, maxy = extracts["bbox"]
    assert minx < maxx
    assert miny < maxy
    slice_bbox = shapely.geometry.box(minx, miny, maxx, maxy)

    admin_level = int(admin_level_slug.replace("admin-level-", ""))
    logging.info(f"Reading level {admin_level} admin boundaries")
    admin_areas: gpd.GeoDataFrame = gpd.read_parquet(admin_areas_path)
    area_unique_id_col = f"GID_{admin_level}"
    area_name_col = f"NAME_{admin_level}"
    # retain the names of larger, encompassing administrative units
    contextual_name_cols = [f"NAME_{i}" for i in range(0, admin_level)]
    admin_areas = admin_areas[[area_name_col, area_unique_id_col, *contextual_name_cols, "geometry"]]
    slice_admin_areas = admin_areas[admin_areas.intersects(slice_bbox)]
    logging.info(f"Found {len(slice_admin_areas)} admin areas intersecting slice bbox")

    if data_to_aggregate.empty:
        logging.info("No damage data; writing admin areas contained slice but no damages")
        slice_admin_areas.to_parquet(aggregated_path)
        sys.exit(0)

    logging.info("Spatially joining input to admin areas")
    # specify left to keep all admin areas for slice (even if they do not appear in the hazard data)
    joined = gpd.sjoin(slice_admin_areas, data_to_aggregate, how="left")
    logging.info(f"{joined.shape=}")

    logging.info(f"Aggregating to {len(slice_admin_areas)} admin boundaries")
    # we only want to groupby certain columns... use the column_regex for this
    columns_to_aggregate = [col for col in joined.columns if re.match(column_regex, col)]
    grouped = joined.groupby(by=area_unique_id_col)[columns_to_aggregate]
    agg_func: str = agg_func_slug.replace("agg-", "")
    aggregated = getattr(grouped, agg_func)()
    aggregated_with_geometry = slice_admin_areas.merge(aggregated, on=area_unique_id_col)

    logging.info(f"Writing {aggregated_with_geometry.shape} to disk")
    aggregated_with_geometry.to_parquet(aggregated_path)

    logging.info("Done")
