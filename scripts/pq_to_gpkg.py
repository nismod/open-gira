#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
import sys

import geopandas as gpd
import pyogrio


GPKG_EXT = "gpkg"
PARQUET_EXT = ("pq", "parq", "parquet")
GEOPARQUET_EXT = ("gpq", "geoparq", "geoparquet")


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
    gpkg_path = f"{stem}.{GPKG_EXT}"

    # vectorised file io with pyogrio, ~50x faster than gpd.to_file via fiona
    pyogrio.write_dataframe(gdf, gpkg_path)

    return gpkg_path


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    if len(sys.argv) < 2:
        sys.exit(
            f"usage:\npython {sys.argv[0]} <geoparquet_file_path_1> <geoparquet_file_path_2> ...\n"
            "a geopackage file will be written beside each input parquet file"
        )
    parquet_files: list = sys.argv[1:]

    for infile in parquet_files:
        file_ext = infile.lower().split(".")[-1]
        if not (file_ext in PARQUET_EXT or file_ext in GEOPARQUET_EXT):
            raise ValueError("{infile=} does not appear to be geoparquet")
        gpkg_path = geoparquet_to_geopackage(infile)
        logging.info(gpkg_path)
