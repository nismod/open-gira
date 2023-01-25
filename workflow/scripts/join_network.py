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
    try:
        node_slice_files = snakemake.input["nodes"]  # type: ignore
        edge_slice_files = snakemake.input["edges"]  # type: ignore
        nodes_output_file = snakemake.output["nodes"]  # type: ignore
        edges_output_file = snakemake.output["edges"]  # type: ignore
    except NameError:
        # not sure of an elegant way to handle the two lists of input filenames
        sys.exit("please invoke via snakemake")

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    logging.info("Joining network node slices")
    nodes = concat_slices(node_slice_files)

    logging.info("Joining network edge slices")
    edges = concat_slices(edge_slice_files)

    # drop the ids we used on a per-slice basis
    nodes = nodes.drop(["node_id"], axis="columns")
    edges = edges.drop(["edge_id", "from_node_id", "to_node_id"], axis="columns")

    network = snkit.network.Network(edges=edges, nodes=nodes)

    # relabel with network-wide ids prior to adding topology
#    logging.info("Labelling edges and nodes with ids")
#    network = snkit.network.add_ids(network)

#    logging.info("Labelling edge ends with node ids")
#    network = snkit.network.add_topology(network)

    # TODO: adding the topology and component ids needs accelerating, 9M rows takes over a day

    # TODO: check what takes the time, is it identifying the components or labelling them?

    # perhaps, do it for each slice, and then at this step, check which slice components
    # join neighbouring slice components and relabel components as such
    # N.B. this will require using the {start|end}_node_reference=NaN nodes
    # N.B. these should only be at bbox edges, but appear to be all over slices

    # another option is to use igraph rather than networkx in snkit to identify the components

#    logging.info("Labelling edges and nodes with network component ids")
#    network = snkit.network.add_component_ids(network)

    logging.info("Writing network to disk")
    network.nodes.to_parquet(nodes_output_file)
    network.edges.to_parquet(edges_output_file)
