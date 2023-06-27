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

from open_gira.io import concat_geoparquet
from open_gira.utils import natural_sort


if __name__ == "__main__":
    try:
        slice_files = snakemake.input  # type: ignore
        output_file = snakemake.output[0]  # type: ignore
    except NameError:
        slice_files = sys.argv[1:-1]
        output_file = sys.argv[-1]

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    logging.info(f"Reading {len(slice_files)=} files")

    joined: gpd.GeoDataFrame = concat_geoparquet(slice_files)

    folder_path = os.path.dirname(os.path.abspath(output_file))
    if not os.path.exists(folder_path):
        os.path.makedirs(folder_path)

    logging.info(f"Writing {joined.shape=} to {output_file}")

    joined.to_parquet(output_file)
