"""Takes a region and finds the box_ids within it, exports as a list.

# TODO delete this script if region-boxes.smk is unnecessary

"""
import json
import os

import geopandas
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.polygon import Polygon


if __name__ == "__main__":
    global_boxes_metadata_path = snakemake.input["global_boxes_metadata"]  # type: ignore
    global_boxes_path = snakemake.input["global_boxes"]  # type: ignore
    basin_geometry_path = snakemake.input["basin_geometry"]  # type: ignore
    region_str = snakemake.wildcards.REGION  # type: ignore
    output_dir = snakemake.config["output_dir"]  # type: ignore
    box_ids_path = snakemake.output["boxes_in_region"]  # type: ignore
    box_ids_with_assets_path = snakemake.output["boxes_in_region_with_assets"]  # type: ignore

    regions: geopandas.GeoDataFrame = geopandas.read_file(basin_geometry_path)

    # check the regions use -180 < longitude <= 180
    assert max(regions.geometry.bounds.maxx) <= 180

    # select a region
    region_polygon = regions[regions.name==region_str].geometry.squeeze()
    if not isinstance(region_polygon, (Polygon, MultiPolygon)):
        raise ValueError(f"could not extract suitable region geometry: {region_polygon=} {regions=} {region_str=}")

    # read in the global grid
    global_grid: geopandas.GeoDataFrame = geopandas.read_file(global_boxes_path)
    # check the grid uses -180 < longitude <= 180
    assert max(global_grid.geometry.bounds.maxx) <= 180

    # find box_ids of boxes intersecting region polygon
    boxes_in_region: list[str] = global_grid[global_grid.intersects(region_polygon)].box_id.to_list()

    # write out all boxes within this region
    with open(box_ids_path, "w") as fp:
        json.dump(boxes_in_region, fp)

    # create a second list of box ids -- those containing part of a country,
    # with some infrastructure inside them
    boxes_in_region_with_assets = []

    with open(global_boxes_metadata_path, "r") as fp:
        box_country_dict = json.load(fp)["box_country_dict"]

    for box_id in boxes_in_region:
        box = box_id.replace("box_", "")
        # TODO: ideally, we'd take this filename logic and the box_id loop to
        # rule level -- that way, each box could be a different job and we'd
        # be closer to only processing boxes we're interested in
        # for now, leave as-is to simplify creating the "boxes in the region
        # with infrastructure" list
        gridfinder_box_path = os.path.join(
            output_dir,
            "processed",
            "power",
            box,
            f"nodes_{box}.parquet",
        )
        if os.path.exists(gridfinder_box_path):
            gridfinder_box = geopandas.read_parquet(gridfinder_box_path)
            if len(gridfinder_box) > 0:
                # infrastructure in box
                boxes_in_region_with_assets.append(box_id)

    with open(box_ids_with_assets_path, "w") as fp:
        json.dump(boxes_in_region_with_assets, fp)
