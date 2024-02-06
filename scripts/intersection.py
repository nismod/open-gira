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

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import pyproj
from pyproj import Geod
from snail.core.intersections import split_linestring
from snail.core.intersections import get_cell_indices as get_cell_indicies_of_midpoint
from tqdm import tqdm

from open_gira import fields


def associate_raster(df, fname, band_number=1) -> pd.Series:
    """
    For each split geometry, lookup the relevant raster value. Cell indicies
    must have been previously calculated and stored as fields.RASTER_{I,J}.

    N.B. This will pass through no data values from the raster (no filtering).

    Args:
        df (pd.DataFrame): Table of features, each with cell indicies pertaining
            to relevant raster pixel. Indicies must be stored under columns with
            names referenced by fields.RASTER_I and fields.RASTER_J
        fname (str): Filename of raster file to read data from
        band_number (int): Which band of the raster file to read

    Returns:
        pd.Series: Series of raster values, with same row indexing as df.
    """

    with rasterio.open(fname) as dataset:

        band_data: np.ndarray = dataset.read(band_number)

        # 2D numpy indexing is j, i (i.e. row, column)
        return pd.Series(index=df.index, data=band_data[df[fields.RASTER_J], df[fields.RASTER_I]])


def write_empty_files(columns, outputs_path):
    try:
        empty_geodf = gpd.GeoDataFrame(
            columns=columns,
            geometry="geometry",
            crs=pyproj.CRS.from_user_input(4326)
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
    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)
    tqdm.pandas()
    try:
        network_edges_path: str = snakemake.input.network
        raster_paths: list[str] = snakemake.input.tif_paths
        copy_raster_values: bool = snakemake.params.copy_raster_values
        output_path: str = snakemake.output.geoparquet
    except NameError:
        sys.exit("Please run from snakemake")

    if not isinstance(raster_paths, list):
        raise ValueError(f"input tif_paths object is not a list, quitting.")

    if len(raster_paths) == 0:
        raise ValueError("The list of TIFF files is empty, quitting.")

    raster_basenames = [
        re.sub("\\.tif$", "", os.path.basename(tif)) for tif in raster_paths
    ]

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
    logging.info(f"{raster_width=} {raster_height=}")

    if len(raster_paths) > 1:
        # Check all raster files use the same grid
        logging.info("Checking raster grid consistency")
        for raster_path in tqdm(raster_paths[1:]):
            with rasterio.open(raster_path) as raster:
                if (
                    raster_width != raster.width
                    or raster_height != raster.height
                    or raster_transform != list(raster.transform)
                ):
                    raise AttributeError(
                        (
                            f"Raster attribute mismatch in file {raster_path}:\n"
                            f"Height: expected={raster_height}; actual={raster.height}\n"
                            f"Width: expected={raster_width}; actual={raster.width}\n"
                            f"Transform equal? {'True' if list(raster.transform) == raster_transform else 'False'}"
                        )
                    )

    # Split edges
    logging.info("Split edges")
    core_splits = []
    for i in tqdm(range(len(core_edges))):
        # split edge
        splits = split_linestring(
            core_edges.geometry[i], raster_width, raster_height, raster_transform
        )
        # add to collection
        for s in splits:
            # splitting sometimes returns zero-length linestrings on edge of raster
            # see below for example linestring on eastern (lon=70W) extent of box
            # (Pdb) geometry.coords.xy
            # (array('d', [-70.0, -70.0]), array('d', [18.445832920952196, 18.445832920952196]))
            # this split geometry has: j = raster_width
            # however j should be in range: 0 <= j < raster_width
            # as a hacky workaround, drop any splits with length 0
            # do we need a nudge off a cell boundary somewhere when performing the splits in snail?
            if not s.length == 0:
                new_row = core_edges.iloc[i].copy()
                new_row.geometry = s
                core_splits.append(new_row)

    core_splits = gpd.GeoDataFrame(core_splits)

    logging.info("Split %d edges into %d pieces", len(core_edges), len(core_splits))

    logging.info("Find indices")

    def cell_indicies_of_split_geometry(geometry, *args, **kwargs) -> pd.Series:
        """
        Given a geometry, find the cell index (i, j) of its midpoint for the
        enclosing raster parameters.

        N.B. There is no checking whether a geometry spans more than one cell.
        """

        # integer indicies
        i, j = get_cell_indicies_of_midpoint(geometry, raster_height, raster_width, raster_transform)

        # die if we're out of bounds somehow
        assert 0 <= i < raster_width
        assert 0 <= j < raster_height

        # return a series with labels so we can unpack neatly into two dataframe columns
        return pd.Series(index=(fields.RASTER_I, fields.RASTER_J), data=[i, j])

    core_splits = pd.concat(
        [
            core_splits,
            core_splits.geometry.progress_apply(
                cell_indicies_of_split_geometry,
                result_type='expand',
            )
        ],
        axis="columns"
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
            raster_data[f"{fields.HAZARD_PREFIX}{raster_basenames[i]}"] = associate_raster(core_splits, raster_paths[i])

        raster_data = pd.DataFrame(raster_data)
        core_splits = pd.concat([core_splits, raster_data], axis="columns")
        assert len(raster_data) == len(core_splits)

    logging.info(f"Write data {core_splits.shape=}")
    core_splits.to_parquet(output_path)

    logging.info("Write data without geometry")
    pd.DataFrame(core_splits.drop(columns=["geometry"])).to_parquet(
        re.sub("\\.geoparquet$", ".parquet", output_path, re.IGNORECASE)
    )

    logging.info("Done.")
