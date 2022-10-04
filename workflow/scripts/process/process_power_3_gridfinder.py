"""Indexes all gridfinder values"""
import json
import os
import sys
import time
import warnings

import fiona
import geopandas as gpd
import numpy as np
from shapely.geometry import shape
from tqdm import tqdm

warnings.filterwarnings("ignore", category=DeprecationWarning)
from .process_power_functions import idxbox

try:
    output_dir = snakemake.params["output_dir"]  # type: ignore
except:
    output_dir = sys.argv[1]


if __name__ == "__main__":

    # preliminary function variables
    with open(
        os.path.join(output_dir, "power_processed", "world_boxes_metadata.json"), "r"
    ) as filejson:
        world_boxes_metadata = json.load(filejson)
    boxlen = world_boxes_metadata["boxlen"]
    lat_max = world_boxes_metadata["lat_max"]
    lon_min = world_boxes_metadata["lon_min"]
    num_cols = world_boxes_metadata["num_cols"]
    tot_boxes = world_boxes_metadata["tot_boxes"]

    s = time.time()

    features = []
    with fiona.open(
        os.path.join(output_dir, "input", "gridfinder", "grid.gpkg")
    ) as src:

        for jj, feature in tqdm(
            enumerate(src), desc="loading grid.gpkg features", total=len(src)
        ):
            # gridfinder GeoPackage stores an "fid" which GeoPandas ignores
            # and fiona reads as "id", not to feature['properties']
            # see https://github.com/geopandas/geopandas/issues/1035
            geom = shape(feature["geometry"])

            features.append(
                {
                    "source_id": feature["id"],
                    "source": feature["properties"]["source"],
                    "geometry": geom,
                }
            )

    gdf = gpd.GeoDataFrame(features)
    print("time for grid.gpkg processing: ", round((time.time() - s) / 60, 2), " mins")

    centres = gdf.centroid
    print(f"centres found: {round((time.time() - s)/60, 2)}. Finding indices...")
    all_lon, all_lat = np.array(centres.x), np.array(centres.y)

    print(f"indices found: {round((time.time() - s)/60, 2)}. Saving...")
    gdf["box_id"] = idxbox(
        all_lat, all_lon, boxlen, lat_max, num_cols, lon_min, tot_boxes
    )

    for box_id, gdf_box in tqdm(
        gdf.groupby("box_id"),
        desc="saving gridfinder gpkg",
        total=len(gdf["box_id"].unique()),
    ):
        all_boxes_path = os.path.join(
            output_dir, "power_processed", "all_boxes", f"{box_id}"
        )
        p = os.path.join(all_boxes_path, f"gridfinder_{box_id}.gpkg")
        gdf_box.to_file(p, driver="GPKG")

    with open(
        os.path.join(output_dir, "power_processed", "world_boxes_metadata.json"), "r"
    ) as filejson:
        tot_boxes = json.load(filejson)["tot_boxes"]

    cols = ["source_id", "source", "box_id", "geometry"]
    for id in tqdm(
        range(int(float(tot_boxes))), desc="empty folders", total=float(tot_boxes)
    ):  # create empty ones
        all_boxes_g_file = os.path.join(
            output_dir,
            "power_processed",
            "all_boxes",
            f"box_{id}",
            f"gridfinder_box_{id}.gpkg",
        )
        if not os.path.exists(all_boxes_g_file):
            empty_file = gpd.GeoDataFrame(columns=cols)
            empty_file.loc[0, :] = [None] * len(cols)
            empty_file["box_id"] = f"box_{id}"
            empty_file.to_file(all_boxes_g_file, driver="GPKG")

    print("gpkg saving: ", round((time.time() - s) / 60, 2), " mins")
