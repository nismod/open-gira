import logging
import sys
import warnings
from typing import Iterable

import geopandas as gpd
import snkit

from open_gira.io import concat_geoparquet
from open_gira.utils import natural_sort


def concat_slices(slice_files: Iterable[str]) -> gpd.GeoDataFrame:
    concatenated: gpd.GeoDataFrame = concat_geoparquet(slice_files)
    logging.info(f"{len(concatenated)=}")

    # drop_duplicates on a GeoDataFrame can be extremely slow
    # instead use something easily hashable
    wkt: pd.Series = concatenated["geometry"].apply(lambda geom: geom.wkt)
    deduplicated = concatenated.loc[wkt.drop_duplicates().index, :]
    deduplicated = deduplicated.reset_index(drop=True)
    logging.info(f"{len(deduplicated)=}")

    return deduplicated


if __name__ == "__main__":
    node_slice_files = snakemake.input["nodes"]
    edge_slice_files = snakemake.input["edges"]
    nodes_output_file = snakemake.output["nodes"]
    edges_output_file = snakemake.output["edges"]

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    logging.info("Joining network node slices")
    nodes = concat_slices(node_slice_files)

    logging.info("Joining network edge slices")
    edges = concat_slices(edge_slice_files)

    network = snkit.network.Network(edges=edges, nodes=nodes)

    # ensure we haven't reused any ids
    assert len(network.nodes) == len(network.nodes.id.unique())
    assert len(network.edges) == len(network.edges.id.unique())

    logging.info("Labelling edge ends with from/to node ids")
    network = snkit.network.add_topology(network)

    logging.info("Labelling edges and nodes with network component ids")
    network = snkit.network.add_component_ids(network)

    logging.info("Writing network to disk")
    network.nodes.to_parquet(nodes_output_file)
    network.edges.to_parquet(edges_output_file)
