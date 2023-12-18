"""
Create a network from plants, targets and gridfinder line data
"""

import logging
import os
from typing import Callable
import sys

import geopandas as gpd
import pandas as pd
from pyproj import Geod
from shapely.geometry import Point, LineString

from open_gira.grid import weighted_allocation


MAX_LINK_DISTANCE = 20_000


def node_to_edge_distance(node: Point, edge: LineString, geoid: Geod = Geod(ellps="WGS84")) -> float:
    """
    Calculate the great circle distance from node to nearest point on edge,
    given a geoid.
    """

    # find the location, p, on the linestring nearest the node
    p: Point = edge.interpolate(edge.project(node))

    _, _, distance = geoid.inv(node.x, node.y, p.x, p.y)

    return distance


def get_proximity_condition(asset_type: str, edge_limit: float | int) -> Callable:
    """
    Return a closure (enclosed function) to check if some node and asset
    pair satisfy a proximity condition.

    Args:
        asset_type: Assets to check
        edge_limit: Maximum distance from node to edge

    Returns:
        Closure function
    """

    geoid = Geod(ellps="WGS84")

    def closure(node: pd.Series, edge: pd.Series) -> bool:
        """Return true if node is a target and within distance limit."""
        return (node.asset_type == asset_type) and (node_to_edge_distance(node.geometry, edge.geometry, geoid) < edge_limit)

    return closure


if __name__ == "__main__":
    plants_path: str = snakemake.input.plants
    targets_path: str = snakemake.input.targets
    gridfinder_path: str = snakemake.input.gridfinder
    snkit_processes: int = snakemake.threads
    nodes_path: str = snakemake.output.nodes
    edges_path: str = snakemake.output.edges
    grid_hull_path: str = snakemake.output.grid_hull

    # set the environment variable governing how many threads snkit will use for parallel operations
    os.environ["SNKIT_PROCESSES"] = str(snkit_processes)
    # delayed import, must wait until environment variable is present
    import snkit
    import snkit.network

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    logging.info("Read node data (plants and targets)")
    plants = gpd.read_parquet(plants_path)
    plants = plants[["source_id", "power_mw", "asset_type", "geometry"]]

    targets = gpd.read_parquet(targets_path)
    targets.geometry = targets.geometry.centroid
    # should not actually reproject, but CRS metadata must match exactly for concat
    targets = targets.to_crs(plants.crs)

    targets["target_id"] = targets["id"]

    nodes = gpd.GeoDataFrame(
        pd.concat([plants, targets], ignore_index=True),
        crs=plants.crs
    )

    logging.info("Read edge data (transmission and distribution lines)")
    edges = gpd.read_parquet(gridfinder_path).reset_index(names="source_id")
    edges["asset_type"] = "transmission"
    edges = edges[["source_id", "asset_type", "geometry", "source"]]

    network = snkit.network.Network(nodes, edges)
    logging.info(f"Raw network contains: {len(nodes)} nodes and {len(edges)} edges")

    if len(nodes) == 0 or len(edges) == 0:
        logging.info("No viable network, writing empty files to disk")
        logging.info("Writing network to disk")
        empty_gdf = gpd.GeoDataFrame({"geometry": []}, crs=4326)
        empty_gdf.to_parquet(edges_path)
        empty_gdf.to_parquet(nodes_path)

        logging.info("Writing network's convex hull to disk")
        empty_gdf.to_file(grid_hull_path)

        sys.exit(0)

    logging.info("Cleaning network")

    logging.info("Split multilinestrings")
    network = snkit.network.split_multilinestrings(network)

    # Connect targets
    logging.info("Connect targets to grid")
    network.nodes["id"] = range(len(network.nodes))
    network = snkit.network.link_nodes_to_nearest_edge(
        network,
        get_proximity_condition("target", MAX_LINK_DISTANCE)
    )
    network.nodes.loc[network.nodes.id.isnull(), "asset_type"] = "conn_target"

    # Connect power plants
    logging.info("Connect power plants to grid")
    network.nodes["id"] = range(len(network.nodes))
    network = snkit.network.link_nodes_to_nearest_edge(
        network,
        get_proximity_condition("source", MAX_LINK_DISTANCE)
    )
    network.nodes.loc[network.nodes.id.isnull(), "asset_type"] = "conn_source"

    # Add nodes at line endpoints
    logging.info("Add nodes at edge ends")
    network.nodes["id"] = range(len(network.nodes))
    network = snkit.network.add_endpoints(network)
    network.nodes.loc[network.nodes.id.isnull(), "asset_type"] = "intermediate"

    logging.info(f"Cleaned network contains: {len(network.nodes)} nodes and {len(network.edges)} edges")

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

    targets = network.nodes[network.nodes.asset_type == "target"]

    # check we have a GDP figure for every target, even if it's 0
    if any(targets.gdp.isna()):
        raise ValueError("Cannot allocate power without a GDP figure to weight by for every target.")

    # if there's no gdp data available at all, use the population as a weight
    # this weight must then be used when allocating after failure in grid_disruption.py
    if targets.gdp.sum() == 0:
        weight_col="population"
    else:
        weight_col="gdp"

    powered_targets: pd.DataFrame = weighted_allocation(
        network.nodes,
        variable_col="power_mw",
        weight_col=weight_col,
        component_col="component_id",
        asset_col="asset_type",
        source_name="source",
        sink_name="target",
    )

    # merge the target power allocation back into the nodes table
    # first get the target power in a table of the same length as network.nodes
    intermediate = network.nodes[["id", "power_mw"]].merge(
        powered_targets[["id", "power_mw"]],
        how="left",
        on="id",
        suffixes=["_x", ""]
    )
    # then use combine_first to overwrite the NaN target power values in network.nodes
    network.nodes = network.nodes.combine_first(intermediate[["power_mw"]])

    # check power has been allocated appropriately
    for c_id, c_nodes in network.nodes.groupby(network.nodes.component_id):

        power_balance: float = c_nodes.power_mw.sum()
        if power_balance > 1E-6:

            # in the case where there is a generator in a component
            # but no targets, the power balance will be positive
            target_mask: pd.Series = c_nodes.asset_type == "target"
            n_targets: int = len(c_nodes[target_mask])
            if n_targets != 0:

                # sometimes targets don't have population, in which case they won't be allocated power
                population: float = c_nodes[target_mask].population.sum()
                if population != 0:

                    raise ValueError(
                        f"Network component {c_id} has {n_targets} targets, \n"
                        f"positive {population} population but a non-zero \n"
                        f"({power_balance} MW) power balance"
                    )

    logging.info("Writing network to disk")
    network.edges.to_parquet(edges_path)
    network.nodes.to_parquet(nodes_path)

    logging.info("Writing network's convex hull to disk")
    grid_hull = gpd.GeoDataFrame({"geometry": [network.edges.geometry.unary_union.convex_hull]})
    grid_hull.to_file(grid_hull_path)
