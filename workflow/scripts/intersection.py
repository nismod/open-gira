#!/usr/bin/env python
# coding: utf-8
"""Split network edges along grid cells, join hazard values.
"""
import glob
import logging
import os
import re
import sys

import geopandas
import pandas
import rasterio
from pyproj import Geod
from snail.core.intersections import get_cell_indices, split_linestring
from tqdm import tqdm


def main(network_edges_path, hazard_tifs, output_path):
    """
    Split the entries in network_edges_path according to the cells they occupy in the
    grids of the hazard_dir. Write the results to a .geoparquet file (output_path)
    and a similarly-named .parquet file.

    Parameters
    ----------
    network_edges_path (str): Path to a .geoparquet file with network data
    attrs (str|List[str]): attribute/s to copy from the original rows when split
    hazard_dir (List[str]): list of hazard raster files whose values should be combined with the network
    output_path (str): .geoparquet version of the path to write to. A .parquet version is also written

    Returns
    -------
    (void)
    """
    hazard_tifs_basenames = [
        re.sub("\\.tif$", "", os.path.basename(tif)) for tif in hazard_tifs
    ]

    # Read metadata for a single raster
    logging.info("Determining raster grid properties")
    with rasterio.open(hazard_tifs[0]) as dataset:
        raster_width = dataset.width
        raster_height = dataset.height
        raster_transform = list(dataset.transform)

    # Check all raster files use the same grid
    logging.info("Checking raster grid consistency")
    for hazard in tqdm(hazard_tifs[1:]):
        with rasterio.open(hazard) as raster:
            if (
                raster_width != raster.width
                or raster_height != raster.height
                or raster_transform != list(raster.transform)
            ):
                raise AttributeError(
                    (
                        f"Raster attribute mismatch in file {hazard}:\n"
                        f"Height: expected={raster_height}; actual={raster.height}\n"
                        f"Width: expected={raster_width}; actual={raster.width}\n"
                        f"Transform equal? {'True' if list(raster.transform) == raster_transform else 'False'}"
                    )
                )

    # Read network edges
    logging.info("Read edges")
    try:
        core_edges = geopandas.read_parquet(network_edges_path)
    except ValueError:
        logging.info("No data in geometry column, writing empty files.")
        columns = [
            "geometry",
            *pandas.read_parquet(network_edges_path).columns,
            "length_km",
            *hazard_tifs_basenames,
        ]
        write_empty_files(columns, output_path)
        return

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
            new_row = core_edges.iloc[i].copy()
            new_row.geometry = s
            core_splits.append(new_row)

    core_splits = geopandas.GeoDataFrame(core_splits)

    logging.info("Split %d edges into %d pieces", len(core_edges), len(core_splits))

    logging.info("Find indices")

    def get_indices(geom):
        x, y = get_cell_indices(geom, raster_width, raster_height, raster_transform)
        x = x % raster_width
        y = y % raster_height
        return [x, y]

    core_splits["cell_index"] = core_splits.geometry.progress_apply(get_indices)

    logging.info("Segment length")
    geod = Geod(ellps="WGS84")
    core_splits["length_km"] = (
        core_splits.geometry.progress_apply(geod.geometry_length) / 1e3
    )

    logging.info("Add hazard values")
    for i in tqdm(range(len(hazard_tifs))):
        associate_raster(core_splits, hazard_tifs_basenames[i], hazard_tifs[i])

    logging.info("Write data")
    core_splits.to_parquet(output_path)

    logging.info("Write data without geometry")
    pandas.DataFrame(core_splits.drop(columns=["geometry"])).to_parquet(
        re.sub("\\.geoparquet$", ".parquet", output_path, re.IGNORECASE)
    )

    logging.info("Done.")


def associate_raster(df, key, fname, band_number=1):
    with rasterio.open(fname) as dataset:
        band_data = dataset.read(band_number)
        df[key] = df.cell_index.apply(lambda i: band_data[i[1], i[0]])


def write_empty_files(columns, outputs_path):
    try:
        empty_geodf = geopandas.GeoDataFrame(columns=columns, geometry="geometry")
    except ValueError:
        raise ValueError("Empty dataframe must contain a geometry column")
    logging.info("Write data")
    empty_geodf.to_parquet(outputs_path)
    logging.info("Write data without geometry")
    pandas.DataFrame(empty_geodf.drop(columns=["geometry"])).to_parquet(
        re.sub("\\.geoparquet$", ".parquet", outputs_path, re.IGNORECASE)
    )


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    tqdm.pandas()
    try:
        network_edges_path = snakemake.input["network"]  # type: ignore
        hazard_dir = snakemake.input["tifs"]  # type: ignore
        output_path = snakemake.output["geoparquet"]  # type: ignore
    except NameError:
        print(sys.argv)
        (network_edges_path, hazard_dir, output_path) = sys.argv[1:]
        # network_edges_path = '../../results/geoparquet/tanzania-mini_filter-highway-core/slice-2.geoparquet'
        # output_path = '../../results/test.geoparquet'
        # hazard_dir = '../../results/input/hazard-aqueduct-river/tanzania-mini'

    tifs = glob.glob(os.path.join(hazard_dir, "*.tif"))

    if len(tifs) == 0:
        raise ValueError(
            (
                f"The list of hazard .tif files is empty. Check they were downloaded to "
                f"{hazard_dir}"
            )
        )

    # print(f"hazard_dir={hazard_dir}")
    # print(f"tifs={tifs}")
    main(
        network_edges_path=network_edges_path, hazard_tifs=tifs, output_path=output_path
    )
