"""
Subset IBTrACS (historic) or STORM (synthetic) storm tracks to a slicing box.
"""

import os

import geopandas as gpd
import shapely


# expand the slice by this many degrees
TRACK_SLICING_BUFFER_DEG = 1


if __name__ == "__main__":

    global_tracks_path = snakemake.input.global_tracks
    grid_bbox_path = snakemake.input.grid_bbox
    sliced_tracks_path = snakemake.output.sliced_tracks
    box_number = snakemake.wildcards.BOX

    with open(grid_bbox_path, "r") as fp:
        bbox_dict = fp.read()
    box_geom = shapely.geometry.polygon.Polygon(bbox_dict)

    tracks = gpd.read_parquet(global_tracks_path)

    # take all storm track points intersecting a buffer of the bounding box
    # adding the buffer reduces the likelihood we construct wind fields where a
    # storm suddenly appears well inside the domain
    tracks = tracks[tracks.intersects(box_geom.buffer(TRACK_SLICING_BUFFER_DEG))]

    os.makedirs(os.path.dirname(sliced_tracks_path), exist_ok=True)
    tracks.to_parquet(sliced_tracks_path)
