"""
Takes a list of geoparquet files containing geodataframes and joins them

    file1.geoparguet
    id obs1 obs2 geometry
    0  a    b    geom0
    1  c    d    geom1

    file2.geoparguet
    id obs1 obs2 geometry
    0  A    B    GEOM0
    1  C    D    GEOM1

    joined.geoparguet
    id obs1 obs2 geometry
    0  a    b    geom0
    1  c    d    geom1
    2  A    B    GEOM0
    3  C    D    GEOM1

Usage:
    python join_data.py [FILE] [output]

Example:
    python join_data.py file1.geoparguet file2.geoparguet joined.geoparguet
"""

import logging
import os
import sys
import warnings

import geopandas as gpd
import pandas
from tqdm import tqdm

from transport.utils import NO_GEOM_ERROR_MSG
from network_components import natural_sort


def append_data(base: gpd.GeoDataFrame, slice_files: list[str]) -> gpd.GeoDataFrame:
    """
    Append GeoDataFrames to one another after deserializing from geoparquet

    Args:
        base (gpd.GeoDataFrame): GeoDataFrame to append others to
        slice_files (list[str]): List of geoparquet file paths

    Returns:
        gpd.GeoDataFrame
    """

    dataframes: list[gpd.GeoDataFrame] = []
    for i, slice_path in tqdm(enumerate(slice_files)):

        try:
            gdf = gpd.read_parquet(slice_path)

        except ValueError as error:
            if NO_GEOM_ERROR_MSG in str(error):
                # if the input parquet file does not contain a geometry column,
                # geopandas will raise a ValueError rather than try to procede. we
                # catch that here, but check the error message - to be more
                # specific than catching and suppressing any ValueError

                # use an empty geodataframe to append instead
                gdf = gpd.GeoDataFrame([])

        dataframes.append(gdf)

    logging.info("Joining slice files")

    # pandas concat of iterable containing GeoDataFrames will return a GeoDataFrame
    return pandas.concat([base, *dataframes])


if __name__ == "__main__":
    try:
        slice_files = snakemake.input  # type: ignore
        output_file = snakemake.output[0]  # type: ignore
    except NameError:
        slice_files = sys.argv[1:-1]
        output_file = sys.argv[-1]

    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    logging.info(f"Reading {len(slice_files)=} files")

    # When getting the input files from snakemake, there is no
    # garantee that they will always in the same order. Sort them for
    # consistency. Makes testing easier.
    slice_files = natural_sort(slice_files)

    try:
        base = gpd.read_parquet(slice_files[-1])
    except ValueError as error:
        if NO_GEOM_ERROR_MSG in str(error):
            base = gpd.GeoDataFrame([])
        else:
            raise error

    joined = append_data(base, slice_files).reset_index(drop=True)

    folder_path = os.path.dirname(os.path.abspath(output_file))
    if not os.path.exists(folder_path):
        os.path.makedirs(folder_path)

    logging.info(f"Writing {joined.shape=} to {output_file}")

    joined.to_parquet(output_file)
