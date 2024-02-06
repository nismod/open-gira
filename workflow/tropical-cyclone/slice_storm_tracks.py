"""
Subset storm tracks by an AOI buffer polygon, taking all points in the track from
first arrival to last departure (including voyages outside the AOI).
"""

import os
import json
import logging
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from tqdm import tqdm


logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)


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
    AOI_points = tracks[tracks.intersects(hull.buffer(track_slicing_buffer_deg))]
    logging.info(f"Found {len(AOI_points.track_id.unique())} tracks passing within buffer")

    # some storms impact a country, wander off outside our AOI, then return some days later
    # want time continuity in our tracks (so we don't interpolate wildly, so our eye speed estimates are reasonable)
    logging.info("Indexing tracks' first arrival to last departure (for time continuity)")
    sliced_tracks_by_track: list[pd.DataFrame] = []
    for track_id in tqdm(AOI_points.track_id.unique()):
        df = AOI_points[AOI_points.track_id == track_id]

        # we need at least two track points to calculate the eye velocity and advective winds
        if len(df) > 1:
            arrival, *_, departure = df.index
            sliced = tracks[tracks.track_id == track_id].loc[arrival: departure]
            assert sliced.timestep.is_monotonic_increasing
            assert np.all(np.diff(sliced.timestep) == 1)
            sliced_tracks_by_track.append(sliced)

    if len(sliced_tracks_by_track) > 0:
        sliced_tracks = pd.concat(sliced_tracks_by_track)
        logging.info(
            f"Filtered to {len(sliced_tracks.track_id.unique())} "
            f"valid (n > 1) tracks, with {len(sliced_tracks)} points"
        )
    else:
        sliced_tracks = gpd.GeoDataFrame({"geometry": []}, crs=4326)
        logging.info("No valid tracks found")

    logging.info("Writing tracks subset to disk")
    os.makedirs(os.path.dirname(sliced_tracks_path), exist_ok=True)
    sliced_tracks.to_parquet(sliced_tracks_path)
