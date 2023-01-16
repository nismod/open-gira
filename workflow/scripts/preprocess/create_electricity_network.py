"""
Creates the network from the plants and targets data
"""

import logging
import multiprocessing
import os
from typing import Callable, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
import snkit
import snkit.network
from pyproj import Geod
from shapely.geometry import Point, LineString

from open_gira.grid import weighted_allocation


def link_node(node, network: snkit.network.Network, condition: Callable) -> Tuple[Point, LineString] | Tuple[None, None]:
    """

    """

    # TODO: type hint node in signature
    logging.info(type(node))

    # for each node, check if nearest edge satisfies condition
    edge = snkit.network.nearest_edge(node.geometry, network.edges)
    if condition is not None and not condition(node, edge):
        return None, None

    # add nodes at points-nearest
    point = snkit.network.nearest_point_on_line(node.geometry, edge.geometry)

    if point != node.geometry:
        # add edges linking
        line = LineString([node.geometry, point])
        return point, line
    else:
        return None, None


def link_nodes_to_nearest_edge(network: snkit.network.Network, condition: None | Callable = None) -> snkit.network.Network:
    """
    Link `network` nodes to nearest edge, subject to `condition`. Parallel
    version of snkit.network.link_nodes_to_nearest_edge

    Args:
        network: Network to operate on
        condition: An optional function which takes arguments of node and edge
            and evaluates to a boolean. If true, link node.

    Returns:
        Network with any new features after linking.
    """

    # build input argument list
    chunk_size: int = np.ceil(len(network.nodes) / os.cpu_count()).astype(int)
    args = [
        (network.nodes.iloc[i: i + chunk_size, :].copy(), network, condition)
        for i in range(0, len(network.nodes), chunk_size)
    ]

    new_features = []
    with multiprocessing.Pool() as pool:
        new_features = pool.starmap(link_node, args)

    # unpack generated features
    new_node_geoms, new_edge_geoms = zip(*new_features)

    new_nodes = snkit.network.matching_gdf_from_geoms(network.nodes, filter(lambda x: x is not None, new_node_geoms))
    all_nodes = snkit.network.concat_dedup([network.nodes, new_nodes])

    new_edges = snkit.network.matching_gdf_from_geoms(network.edges, filter(lambda x: x is not None, new_edge_geoms))
    all_edges = snkit.network.concat_dedup([network.edges, new_edges])

    # split edges as necessary after new node creation
    unsplit = snkit.network.Network(nodes=all_nodes, edges=all_edges)

    return snkit.network.split_edges_at_nodes(unsplit)


if __name__ == "__main__":
    output_dir: str = snakemake.config["output_dir"]
    parallel: bool = snakemake.config["parallelise_by_storm"]
    plants_path: str = snakemake.input.plants
    targets_path: str = snakemake.input.targets
    gridfinder_path: str = snakemake.input.gridfinder
    nodes_path: str = snakemake.output.nodes
    edges_path: str = snakemake.output.edges

    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    logging.info("Read node data (plants and targets)")
    plants = gpd.read_parquet(plants_path)
    plants = plants[["source_id", "power_mw", "asset_type", "geometry"]]

    targets = gpd.read_parquet(targets_path)
    target_nodes = targets[["asset_type", "geometry"]].copy()
    target_nodes.geometry = target_nodes.geometry.centroid
    # should not actually reproject, but CRS metadata must match exactly for concat
    target_nodes = target_nodes.to_crs(plants.crs)

    nodes = gpd.GeoDataFrame(
        pd.concat([plants, target_nodes], ignore_index=True),
        crs=plants.crs
    )

    logging.info("Read edge data (transmission and distribution lines)")
    edges = gpd.read_file(gridfinder_path).reset_index(names="source_id")
    edges["asset_type"] = "transmission"
    edges = edges[["source_id", "asset_type", "geometry", "source"]]

    network = snkit.network.Network(nodes, edges)
    logging.info(f"Raw network contains: {len(nodes)} nodes and {len(edges)} edges")

    logging.info("Cleaning network")

    logging.info("Split multilinestrings")
    network = snkit.network.split_multilinestrings(network)

    geod = Geod(ellps="WGS84")
    edge_limit = 20_000  # meters - TODO test limit

    def node_to_edge_distance(point, line, geod):
        b = line.interpolate(line.project(point))
        _, _, distance = geod.inv(point.x, point.y, b.x, b.y)
        return distance

    if parallel:
        link_function = link_nodes_to_nearest_edge
    else:
        link_function = snkit.network.link_nodes_to_nearest_edge

    # Connect targets
    logging.info("Connect targets to grid")
    network.nodes["id"] = range(len(network.nodes))
    network = link_function(
        network,
        lambda node, edge: (
            (node.asset_type == "target")
            and (node_to_edge_distance(node.geometry, edge.geometry, geod) < edge_limit)
        )
    )
    network.nodes.loc[network.nodes.id.isnull(), "asset_type"] = "conn_target"

    # Connect power plants
    logging.info("Connect power plants to grid")
    network.nodes["id"] = range(len(network.nodes))
    network = link_function(
        network,
        lambda node, edge: (
            (node.asset_type == "source")
            and (node_to_edge_distance(node.geometry, edge.geometry, geod) < edge_limit)
        )
    )
    network.nodes.loc[network.nodes.id.isnull(), "asset_type"] = "conn_source"

    # Add nodes at line endpoints
    # TODO: is this necessary?
    logging.info("Add nodes at edge ends")
    network.nodes["id"] = range(len(network.nodes))
    network = snkit.network.add_endpoints(network)
    network.nodes.loc[network.nodes.id.isnull(), "asset_type"] = "intermediate"

    logging.info(f"Cleaned network contains: {len(network.nodes)} nodes and {len(network.edges)} edges")

    # join econometric and demographic data to network nodes dataframe
    network.nodes = pd.merge(network.nodes, targets[["id", "gdp", "population"]], on="id", how="outer")

    # overwrite ids to int type (grid_disruption.py needs integer ids for fast pandas isin ops)
    network.nodes["id"] = range(len(network.nodes))
    network.edges["id"] = range(len(network.edges))

    # add from/to ids
    logging.info("Infer network topology")
    network = snkit.network.add_topology(network, id_col="id")

    # add connected component ids
    logging.info("Detecting connected components")
    if "component_id" not in network.nodes.columns:
        network = snkit.network.add_component_ids(network)

    logging.info("Allocating generating capacity to targets")
    targets: pd.DataFrame = weighted_allocation(
        network.nodes,
        variable_col="power_mw",
        weight_col="gdp",
        component_col="component_id",
        asset_col="asset_type",
        source_name="source",
        sink_name="target",
    )

    # merge the target power allocation back into the nodes table
    # first get the target power in a table of the same length as network.nodes
    target_power = network.nodes[["id", "power_mw"]].merge(targets[["id", "power_mw"]], how="left", on="id", suffixes=["_x", ""])
    # then use combine_first to overwrite the NaN target power values in network.nodes
    network.nodes = network.nodes.combine_first(target_power[["id", "power_mw"]])

    # check power has been allocated appropriately
    assert network.nodes["power_mw"].sum() < 1E-6

    logging.info("Writing network to disk")
    network.edges.to_parquet(edges_path)
    network.nodes.to_parquet(nodes_path)
