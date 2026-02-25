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


logging.basicConfig(
    format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
)


if __name__ == "__main__":
    global_tracks_path: str = snakemake.input.global_tracks  # noqa: F821
    grid_hull_path: str = snakemake.input.grid_hull  # noqa: F821
    track_slicing_buffer_deg: float = snakemake.config[  # noqa: F821
        "max_track_search_radius_deg"
    ]
    sliced_tracks_path: str = snakemake.output.sliced_tracks  # noqa: F821

    # use a hull to reject tracks which would intersect a bounding box
    # but not impact the grid
    with open(grid_hull_path, "r") as fp:
        data = json.load(fp)

    try:
        (shape_dict,) = data["features"]
    except ValueError:
        logging.info("No network, therefore no intersecting storm tracks")
        os.makedirs(os.path.dirname(sliced_tracks_path), exist_ok=True)
        gpd.GeoDataFrame({"geometry": []}, crs=4326).to_parquet(sliced_tracks_path)
        sys.exit(0)

    hull = shapely.geometry.shape(shape_dict["geometry"])

    tracks = gpd.read_parquet(global_tracks_path)

    logging.info("Subsetting tracks to those which intersect with grid (+ buffer)")
    aoi_points = tracks[tracks.intersects(hull.buffer(track_slicing_buffer_deg))]

    logging.info("Slicing tracks' first arrival to last departure")
    if not aoi_points.index.name:
        aoi_points.index.name = "datetime"
    aoi_grouped = (
        aoi_points.reset_index()
        .groupby("track_id")
        .agg({aoi_points.index.name: ["min", "max"]})
    )
    aoi_grouped.columns = aoi_grouped.columns.droplevel(0)
    aoi_grouped.columns = ["arrival_index", "departure_index"]

    # Filter to tracks that have more than one point in AOI
    valid_aoi_tracks = aoi_grouped[
        aoi_grouped["departure_index"] > aoi_grouped["arrival_index"]
    ]

    if len(valid_aoi_tracks) > 0:
        sliced_tracks_list = []
        tracks_grouped = tracks.groupby("track_id")

        for track_id, row in valid_aoi_tracks.iterrows():
            full_track = tracks_grouped.get_group(track_id)
            sliced_track = full_track.loc[row["arrival_index"] : row["departure_index"]]

            try:
                assert sliced_track.timestep.is_monotonic_increasing
            except AssertionError:
                logging.info(f"{track_id} without monotonically increasing timestep")
                continue
            try:
                assert np.all(np.diff(sliced_track.timestep) == 1)
            except AssertionError:
                logging.info(f"{track_id} has gaps, dropping")
                continue

            sliced_tracks_list.append(sliced_track)

        sliced_tracks = pd.concat(sliced_tracks_list)

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
