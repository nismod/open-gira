# Takes a list of geoparquet files containing geodataframes and joins them

# file1.geoparguet
# id obs1 obs2 geometry
# 0  a    b    geom0
# 1  c    d    geom1

# file2.geoparguet
# id obs1 obs2 geometry
# 0  A    B    GEOM0
# 1  C    D    GEOM1

# joined.geoparguet
# id obs1 obs2 geometry
# 0  a    b    geom0
# 1  c    d    geom1
# 2  A    B    GEOM0
# 3  C    D    GEOM1

# Usage: python join_data.py [FILE] [output]
# Example: python join_data.py file1.geoparguet file2.geoparguet joined.geoparguet

import logging
import sys
import warnings

import geopandas as gpd
import pandas


def append_data(base, slice_files):
    slice_files.pop()
    if len(slice_files) == 0:
        return base

    try:
        gdf = gpd.read_parquet(slice_files[-1])
    except ValueError as error:
        # if the input parquet file does not contain a geometry column, geopandas
        # will raise a ValueError rather than try to procede
        logging.info(f"{error}\n")

        # snakemake requires that output files exist though, so write empty ones
        gdf = gpd.GeoDataFrame([])

    # glue the dataframes together
    # N.B. there is no geopandas concat, so use pandas and then create a new gdf
    base = gpd.GeoDataFrame(pandas.concat([base, gdf]), crs=base.crs)

    return append_data(base, slice_files)


if __name__ == "__main__":
    try:
        slice_files = snakemake.input
        output_file = snakemake.output[0]
    except NameError:
        slice_files = sys.argv[1:-1]
        output_file = sys.argv[-1]

    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    # When getting the input files from snakemake, there is no
    # garantee that they will always in the same order. Sort them for
    # consistency. Makes testing easier.
    slice_files = sorted(slice_files)
    # We're reading the different files as a stack from the top.  Let's
    # reverse the order of files to keep the first file on top.
    slice_files = slice_files[::-1]

    try:
        base = gpd.read_parquet(slice_files[-1])
    except ValueError:
        # if the input parquet file does not contain a geometry column, geopandas
        # will raise a ValueError rather than try to procede
        logging.info("base input file empty... suppressing geopandas exception")

        base = gpd.GeoDataFrame([])

    base = append_data(base, slice_files)
    base = base.reset_index(drop=True)
    base.to_parquet(output_file)
