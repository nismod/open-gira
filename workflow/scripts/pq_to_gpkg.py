#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
import sys

import geopandas as gpd


GPKG_EXT = ".gpkg"
PARQUET_EXT = ".parquet"
GEOPARQUET_EXT = ".geoparquet"


def geoparquet_to_geopackage(parquet_path) -> str:
    gdf = gpd.read_parquet(parquet_path)
    stem, _ = os.path.splitext(parquet_path)
    gpkg_path = stem + GPKG_EXT
    gdf.to_file(gpkg_path)
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
        if not (
            infile.lower().endswith(PARQUET_EXT)
            or infile.lower().endswith(GEOPARQUET_EXT)
        ):
            raise ValueError("{infile=} does not appear to be geoparquet")
        gpkg_path = geoparquet_to_geopackage(infile)
        logging.info(f"{infile} -> {gpkg_path}")
