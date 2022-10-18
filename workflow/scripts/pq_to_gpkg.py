#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
import sys

import geopandas as gpd
import numpy as np


GPKG_EXT = ".gpkg"
PARQUET_EXT = ".geoparquet"
GPKG_INCOMPATIBLE_CLASSES = (np.ndarray,)


def geoparquet_to_geopackage(parquet_path: str) -> str:
    """
    Rewrite geoparquet file as geopackage in same location (new file extension).

    N.B. If geopackage file contains dtypes that are incompatible with geopackage
    implementation in fiona, we will drop these columns and try again.

    Arguments:
        parquet_path (str): Location of geoparquet file

    Returns:
        str: Location of new geopackage file
    """
    gdf = gpd.read_parquet(parquet_path)
    stem, _ = os.path.splitext(parquet_path)
    gpkg_path = stem + GPKG_EXT

    try:
        gdf.to_file(gpkg_path)
    except ValueError:
        # fiona will refuse to serialize some data to geopackage, try dropping

        # hunt for bad dtypes -- and hope first row is representative!
        dtypes = gdf.iloc[0, :].apply(type).values
        bad_indicies = []
        for klass in GPKG_INCOMPATIBLE_CLASSES:
            for i in range(len(dtypes)):
                if dtypes[i] == klass:
                    bad_indicies.append(i)

        logging.warning(f"Dropping {gdf.columns[bad_indicies]} as not compatible with GPKG")
        gdf.drop(columns=gdf.columns[bad_indicies]).to_file(gpkg_path)
    return gpkg_path


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    if len(sys.argv) < 2:
        sys.exit(
            f"usage:\npython {sys.argv[0]} <geoparquet_file_path_1> <geoparquet_file_path_2> ...\n"
            "a geopackage file will be written beside each input parquet file"
        )
    parquet_files: list = sys.argv[1:]

    for infile in parquet_files:
        if not infile.lower().endswith(PARQUET_EXT):
            raise ValueError("{infile=} does not appear to be geoparquet")
        gpkg_path = geoparquet_to_geopackage(infile)
        logging.info(f"{infile} -> {gpkg_path}")
