"""
Subset IBTrACS (historic) or STORM (synthetic) storm tracks by a buffer polygon
"""

import os
import json

import geopandas as gpd
import shapely


# expand the slice by this many degrees
TRACK_SLICING_BUFFER_DEG = 2


if __name__ == "__main__":

    global_tracks_path = snakemake.input.global_tracks
    grid_hull_path = snakemake.input.grid_hull
    sliced_tracks_path = snakemake.output.sliced_tracks

    # use a hull to reject tracks which would intersect a bounding box
    # but not impact the grid
    with open(grid_hull_path, "r") as fp:
        data = json.load(fp)
    shape_dict, = data["features"]
    hull = shapely.geometry.shape(shape_dict["geometry"])

    tracks = gpd.read_parquet(global_tracks_path)

    # adding the buffer reduces the likelihood we construct wind fields where a
    # storm suddenly appears well inside the domain
    tracks = tracks[tracks.intersects(hull.buffer(TRACK_SLICING_BUFFER_DEG))]

    os.makedirs(os.path.dirname(sliced_tracks_path), exist_ok=True)
    tracks.to_parquet(sliced_tracks_path)
