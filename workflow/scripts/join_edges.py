# Takes a list of geoparquet files containing edges and join corresponding
# geodataframes. The geodataframes must have the following columns:
# start_node_longitude start_node_latitude start_node_reference end_node_longitude end_node_latitude end_node_reference geometry

# example join:

# file1.geoparguet
# id obs1 obs2 start_node_longitude start_node_latitude start_node_reference end_node_longitude end_node_latitude end_node_reference geometry
# 0  a    b    0                    0                   0                    1                  1                 1                  geom0
# 1  c    d    1                    1                   1                    2                  1                 2                  geom1

# file2.geoparguet
# id obs1 obs2 start_node_longitude start_node_latitude start_node_reference end_node_longitude end_node_latitude end_node_reference geometry
# 0  A    B    2                    2                   2                    4                  3                 3                  GEOM0
# 1  C    D    4                    3                   3                    4                  4                 4                  GEOM1

# joined.geoparguet
# id obs1 obs2 start_node_longitude start_node_latitude start_node_reference end_node_longitude end_node_latitude end_node_reference geometry
# 0  a    b    0                    0                   0                    1                  1                 1                  geom0
# 1  c    d    1                    1                   1                    2                  1                 2                  geom1
# 0  A    B    2                    2                   2                    4                  3                 3                  GEOM0
# 1  C    D    4                    3                   3                    4                  4                 4                  GEOM1

# Usage: python join_edges [FILE] [output]
# Example: python join_edges file1.geoparguet file2.geoparguet joined.geoparguet

import logging
import sys
import warnings

import geopandas as gpd
import numpy as np
import pandas

from join_data import append_data


def add_custom_node_references(base):
    """
    When converting to .geoparquet we added nodes at the bounding box edges.
    These nodes have no reference. We need to make it easy to identify nodes by
    ensuring that nodes in the same location have the same reference.  We'll
    make it easy on ourselves by giving our inserted nodes negative reference
    numbers.
    """

    # Find start nodes with no reference
    na_start_nodes = (
        base[base.start_node_reference.isna()][
            ["start_node_longitude", "start_node_latitude"]
        ]
        .copy()
        .rename(columns={"start_node_longitude": "lon", "start_node_latitude": "lat"})
    )
    # and end nodes with no reference
    na_end_nodes = (
        base[base.end_node_reference.isna()][
            ["end_node_longitude", "end_node_latitude"]
        ]
        .copy()
        .rename(columns={"end_node_longitude": "lon", "end_node_latitude": "lat"})
    )
    # stitch them together, dropping any duplicate coordinates
    nodes = pandas.concat([na_start_nodes, na_end_nodes]).drop_duplicates()
    # give them ids
    nodes_n = len(nodes)
    nodes["node_reference"] = np.arange(nodes_n)[::-1] - nodes_n

    # merge on against start nodes and fill na values
    base = base.merge(
        nodes,
        left_on=["start_node_longitude", "start_node_latitude"],
        right_on=["lon", "lat"],
        how="left",
    ).drop(columns=["lon", "lat"])
    base.start_node_reference = base.start_node_reference.fillna(base.node_reference)
    base = base.drop(columns="node_reference")

    # merge on against end nodes and fill na values
    base = base.merge(
        nodes,
        left_on=["end_node_longitude", "end_node_latitude"],
        right_on=["lon", "lat"],
        how="left",
    ).drop(columns=["lon", "lat"])
    base.end_node_reference = base.end_node_reference.fillna(base.node_reference)
    base = base.drop(columns="node_reference")

    return base


if __name__ == "__main__":
    try:
        slice_files = snakemake.input  # type: ignore
        output_file = snakemake.output[0]  # type: ignore
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
    base = add_custom_node_references(base)
    base.to_parquet(output_file)
