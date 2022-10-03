"""
Takes a region and finds the box_ids within it, exports as a list.
"""

import json
import os
import sys

import fiona
import geopandas as gpd
from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon


if __name__ == "__main__":

    try:
        global_boxes_metadata_path = snakemake.input["global_boxes_metadata"]
        global_boxes_path = snakemake.input["global_boxes"]
        basin_geometry_path = snakemake.input["basin_geometry"]
        region_str = snakemake.params["region_name"]
        output_dir = snakemake.params["output_dir"]
        box_ids_path = snakemake.output["boxes_in_region"]
        box_ids_with_assets_path = snakemake.output["boxes_in_region_with_assets"]
    except NameError:
        global_boxes_metadata_path = sys.argv[1]
        global_boxes_path = sys.argv[2]
        basin_geometry_path = sys.argv[3]
        region_str = sys.argv[4]
        output_dir = sys.argv[5]
        box_ids_path = sys.argv[6]
        box_ids_with_assets_path = sys.argv[7]

    regions: gpd.GeoDataFrame = gpd.read_file(basin_geometry_path)

    # check the regions use -180 < longitude <= 180
    assert max(regions.geometry.bounds.maxx) <= 180

    # select a region
    region_polygon = regions[regions.name==region_str].geometry.squeeze()
    if not isinstance(region_polygon, (Polygon, MultiPolygon)):
        raise ValueError(f"could not extract suitable region geometry: {region_polygon=} {regions=} {region_str=}")

    # read in the global grid
    global_grid: gpd.GeoDataFrame = gpd.read_file(global_boxes_path)
    # check the grid uses -180 < longitude <= 180
    assert max(global_grid.geometry.bounds.maxx) <= 180

    # find box_ids of boxes intersecting region polygon
    boxes_in_region: list[str] = global_grid[global_grid.intersects(region_polygon)].box_id.to_list()

    # write out all boxes within this region
    with open(box_ids_path, "w") as fp:
        json.dump(boxes_in_region, fp)

    # create a second list of box ids -- those containing part of a country,
    # with some infrastructure inside them
    boxes_in_region_with_assets = boxes_in_region.copy()
    for box_id in boxes_in_region:

        with open(global_boxes_metadata_path, "r") as fp:
            box_country_dict = json.load(fp)["box_country_dict"]

        # empty box will contain None, implying no countries in box
        if not box_country_dict[box_id]:
            boxes_in_region_with_assets.remove(box_id)
            continue

        # TODO: ideally, we'd take this filename logic and the box_id loop to
        # rule level -- that way, each box could be a different job and we'd
        # be closer to only processing boxes we're interested in
        # for now, leave as-is to simplify creating the "boxes in the region
        # with infrastructure" list
        gridfinder_box_path = os.path.join(
            output_dir,
            "power_processed",
            "all_boxes",
            box_id,
            f"gridfinder_{box_id}.gpkg",
        )
        if os.path.exists(gridfinder_box_path):
            gridfinder_box = gpd.read_file(gridfinder_box_path)

            if len(gridfinder_box) == 1 and gridfinder_box["geometry"].iloc[0] == None:
                # no infrastructure in box
                boxes_in_region_with_assets.remove(box_id)
        else:
            # gridfinder file for this box doesn't exist (we haven't produced it)
            boxes_in_region_with_assets.remove(box_id)

    with open(box_ids_with_assets_path, "w") as fp:
        json.dump(boxes_in_region_with_assets, fp)
