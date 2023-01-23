"""
Subset IBTrACS (historic) or STORM (synthetic) storm tracks by a buffer polygon
"""

import os
import json
import logging
import sys

import geopandas as gpd
import shapely


logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)


if __name__ == "__main__":

    global_tracks_path: str = snakemake.input.global_tracks
    grid_hull_path: str = snakemake.input.grid_hull
    track_slicing_buffer_deg: float = snakemake.config["max_track_search_radius_deg"]
    sliced_tracks_path: str = snakemake.output.sliced_tracks

    # use a hull to reject tracks which would intersect a bounding box
    # but not impact the grid
    with open(grid_hull_path, "r") as fp:
        data = json.load(fp)

    try:
        shape_dict, = data["features"]
    except ValueError:
        logging.info("No network, therefore no intersecting storm tracks")
        os.makedirs(os.path.dirname(sliced_tracks_path), exist_ok=True)
        gpd.GeoDataFrame({"geometry": []}, crs=4326).to_parquet(sliced_tracks_path)
        sys.exit(0)

    hull = shapely.geometry.shape(shape_dict["geometry"])

    tracks = gpd.read_parquet(global_tracks_path)

    logging.info("Subsetting tracks to those which intersect with grid (+ buffer)")
    # adding the buffer reduces the likelihood we construct wind fields where a
    # storm suddenly appears well inside the domain
    tracks = tracks[tracks.intersects(hull.buffer(track_slicing_buffer_deg))]

    logging.info("Writing tracks subset to disk")
    os.makedirs(os.path.dirname(sliced_tracks_path), exist_ok=True)
    tracks.to_parquet(sliced_tracks_path)
