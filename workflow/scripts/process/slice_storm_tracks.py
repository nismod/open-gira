"""
Subset IBTrACS (historic) or STORM (synthetic) storm tracks to a slicing box.
"""

import os

import geopandas as gpd


if __name__ == "__main__":

    global_tracks_path = snakemake.input.global_tracks
    global_boxes_path = snakemake.input.global_boxes
    sliced_tracks_path = snakemake.output.sliced_tracks
    box_number = snakemake.wildcards.BOX

    boxes = gpd.read_parquet(global_boxes_path).set_index("box_id")
    box_geom = boxes.loc[f"box_{box_number}", "geometry"]

    tracks = gpd.read_parquet(global_tracks_path)
    tracks = tracks[tracks.intersects(box_geom)]

    os.makedirs(os.path.dirname(sliced_tracks_path), exist_ok=True)
    tracks.to_parquet(sliced_tracks_path)
