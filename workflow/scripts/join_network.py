import logging
import sys
import warnings
from typing import Iterable

import geopandas as gpd
import numpy as np
import pandas

import snkit

from join_data import append_data


def append_slices(slice_files: Iterable[str]) -> gpd.GeoDataFrame:
    # When getting the input files from snakemake, there is no
    # garantee that they will always in the same order. Sort them for
    # consistency. Makes testing easier.

    # We're reading the different files as a stack from the top.  Let's
    # reverse the order of files to keep the first file on top.
    slice_files = sorted(slice_files, reverse=True)

    try:
        base = gpd.read_parquet(slice_files[-1])
    except ValueError:
        # if the input parquet file does not contain a geometry column, geopandas
        # will raise a ValueError rather than try to procede
        logging.info("base input file empty... suppressing geopandas exception")
        base = gpd.GeoDataFrame([])

    concatenated = append_data(base, slice_files)

    deduplicated = concatenated.drop_duplicates("geometry")

    return deduplicated.reset_index(drop=True)


if __name__ == "__main__":
    try:
        node_slice_files = snakemake.input["nodes"]
        edge_slice_files = snakemake.input["edges"]
        nodes_output_file = snakemake.output["nodes"]
        edges_output_file = snakemake.output["edges"]
    except NameError:
        sys.exit("please invoke with snakemake")

    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    logging.info("Joining nodes")
    nodes = append_slices(node_slice_files)

    logging.info("Joining edges")
    edges = append_slices(edge_slice_files)

    network = snkit.network.Network(edges=edges, nodes=nodes)
    logging.info("Labelling edges and nodes with ids")
    network = snkit.network.add_ids(network)

    logging.info("Labelling edge ends with node ids")
    network = snkit.network.add_topology(network)

    logging.info("Labelling edges and nodes with network component ids")
    network = snkit.network.add_component_ids(network)

    logging.info("Writing network to disk")
    network.nodes.to_parquet(nodes_output_file)
    network.edges.to_parquet(edges_output_file)
