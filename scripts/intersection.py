#!/usr/bin/env python
# coding: utf-8
"""
Split network edges along grid cells, identify raster cell for each split
geometry, optionally read hazard values for each split geometry.
"""

import logging
import os
import re
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
import rasterio
import snail.intersection
from pyproj import Geod
from tqdm import tqdm

from open_gira import fields


def write_empty_files(columns, outputs_path):
    try:
        empty_geodf = gpd.GeoDataFrame(
            columns=columns, geometry="geometry", crs=pyproj.CRS.from_user_input(4326)
        )
    except ValueError:
        raise ValueError("Empty dataframe must contain a geometry column")

    logging.info("Write data")
    empty_geodf.to_parquet(outputs_path)
    logging.info("Write data without geometry")
    pd.DataFrame(empty_geodf.drop(columns=["geometry"])).to_parquet(
        re.sub("\\.geoparquet$", ".parquet", outputs_path, re.IGNORECASE)
    )


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )
    tqdm.pandas()
    try:
        network_edges_path: str = snakemake.input.network
        raster_paths: list[str] = snakemake.input.tif_paths
        copy_raster_values: bool = snakemake.params.copy_raster_values
        output_path: str = snakemake.output.geoparquet
    except NameError:
        sys.exit("Please run from snakemake")

    if not isinstance(raster_paths, list):
        # Handle case where input function passes the directory, not the unpacked list of filenames
        raw_folder = Path(raster_paths)
        print(f"{raw_folder=}")
        # where the trimmed tiffs for a given DATASET go
        dataset_folder: Path = raw_folder.parent / snakemake.wildcards.DATASET
        print(f"{dataset_folder=}")

        # file basenames to create (trimmed) full paths for
        raster_paths = []
        for raw_fname in raw_folder.glob("*.tif"):
            raster_paths.append(dataset_folder / raw_fname.name)

    raster_basenames = [
        re.sub("\\.tif$", "", os.path.basename(tif)) for tif in raster_paths
    ]
    if len(raster_paths) == 0:
        raise ValueError("The list of TIFF files is empty, quitting.")

    # Read network edges
    logging.info("Read edges")
    core_edges = gpd.read_parquet(network_edges_path)
    if core_edges.empty:
        logging.info("No data in geometry column, writing empty files.")
        columns = pd.read_parquet(network_edges_path).columns
        write_empty_files(columns, output_path)
        sys.exit(0)

    # Read metadata for a single raster
    logging.info("Determining raster grid properties")
    with rasterio.open(raster_paths[0]) as dataset:
        raster_width = dataset.width
        raster_height = dataset.height
        raster_transform = list(dataset.transform)

    grid = snail.intersection.GridDefinition.from_raster(raster_paths[0])
    logging.info(f"{grid=}")

    if len(raster_paths) > 1:
        # Check all raster files use the same grid
        logging.info("Checking raster grid consistency")
        for raster_path in raster_paths[1:]:
            other_grid = snail.intersection.GridDefinition.from_raster(raster_path)
            if other_grid != grid:
                raise AttributeError(
                    (
                        f"Raster attribute mismatch in file {raster_path}:\n"
                        f"Height: expected={grid.height}; actual={other_grid.height}\n"
                        f"Width: expected={grid.width}; actual={other_grid.width}\n"
                        f"Transform equal? {other_grid.transform == grid.transform}\n"
                        f"Transform expected= {grid.transform}\n"
                        f"Transform actual= {other_grid.transform}\n"
                        f"CRS equal? {other_grid.crs == grid.crs}"
                    )
                )

    # Split edges
    logging.info("Split edges")
    core_splits = snail.intersection.split_linestrings(
        core_edges.reset_index(drop=True), grid
    )

    logging.info("Split %d edges into %d pieces", len(core_edges), len(core_splits))

    logging.info("Find indices")
    core_splits = snail.intersection.apply_indices(
        core_splits, grid, index_i=fields.RASTER_I, index_j=fields.RASTER_J
    )

    logging.info("Calculate split segment lengths")
    geod = Geod(ellps="WGS84")
    core_splits[fields.SPLIT_LENGTH] = (
        core_splits.geometry.progress_apply(geod.geometry_length) / 1e3
    )

    if copy_raster_values:
        # N.B. this loop is the heavy lifting of this script
        # it reads hazard intensity values len(raster_paths) * len(core_splits) times
        logging.info("Adding raster values to split geometries")

        # to prevent a fragmented dataframe (and a memory explosion), add series to a dict
        # and then concat afterwards -- do not append to an existing dataframe
        raster_data: dict[str, pd.Series] = {}

        for i in tqdm(range(len(raster_paths))):
            colname = f"{fields.HAZARD_PREFIX}{raster_basenames[i]}"
            with rasterio.open(raster_paths[i]) as src:
                data = src.read(1)
                raster_data[colname] = snail.intersection.get_raster_values_for_splits(
                    core_splits, data, index_i=fields.RASTER_I, index_j=fields.RASTER_J
                )

        raster_data = pd.DataFrame(raster_data)
        core_splits = pd.concat([core_splits, raster_data], axis="columns")
        assert len(raster_data) == len(core_splits)

    logging.info(f"Write data {core_splits.shape=} {core_splits.columns=}")
    core_splits.to_parquet(output_path)

    logging.info("Write data without geometry")
    pd.DataFrame(core_splits.drop(columns=["geometry"])).to_parquet(
        re.sub("\\.geoparquet$", ".parquet", output_path, re.IGNORECASE)
    )

    logging.info("Done.")
