"""
Creates the network from the plants and targets data
"""

import warnings

import geopandas as gpd
import pandas as pd
import snkit
import snkit.network
from pyproj import Geod
from shapely.errors import ShapelyDeprecationWarning

from open_gira.grid import allocate_power_to_targets


if __name__ == "__main__":
    output_dir = snakemake.config["output_dir"]  # type: ignore
    box_id = snakemake.wildcards.BOX  # type: ignore
    plants_path = snakemake.input.plants  # type: ignore
    targets_path = snakemake.input.targets  # type: ignore
    gridfinder_path = snakemake.input.gridfinder  # type: ignore
    nodes_path = snakemake.output.nodes  # type: ignore
    edges_path = snakemake.output.edges  # type: ignore

    # Nodes
    plants = gpd.read_parquet(plants_path)
    targets = gpd.read_parquet(targets_path)

    plants["id"] = [f"source_{i}_{box_id}" for i in range(len(plants))]
    plants = plants[["id", "source_id", "power_mw", "asset_type", "geometry"]]

    targets["id"] = [f"target_{i}_{box_id}" for i in range(len(targets))]
    target_cols = list(targets.columns)

    target_nodes = targets[["id", "asset_type", "geometry"]].copy()
    target_nodes.geometry = target_nodes.geometry.centroid

    # should not actually reproject, but CRS metadata must match exactly for concat
    target_nodes = target_nodes.to_crs(plants.crs)

    nodes = gpd.GeoDataFrame(
        pd.concat([plants, target_nodes], ignore_index=True),
        crs=plants.crs
    )

    # Edges
    edges = gpd.read_parquet(gridfinder_path)

    edges["asset_type"] = "transmission"

    edges = edges[["source_id", "asset_type", "geometry", "source"]]

    # Process network
    network = snkit.network.Network(nodes, edges)

    if len(edges) > 0:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
            network = snkit.network.split_multilinestrings(network)

            geod = Geod(ellps="WGS84")
            edge_limit = 20_000  # meters - TODO test limit
            def node_to_edge_distance(point, line, geod):
                b = line.interpolate(line.project(point))
                _, _, distance = geod.inv(point.x, point.y, b.x, b.y)
                return distance

            # Connect power plants
            network = snkit.network.link_nodes_to_nearest_edge(
                network,
                lambda node, edge: (
                    (node.asset_type == "source")
                    and (node_to_edge_distance(node.geometry, edge.geometry, geod) < edge_limit)
                )
            )
            network.nodes.loc[network.nodes.id.isnull(), "asset_type"] = "conn_source"

            network.nodes["id"] = network.nodes.reset_index().apply(
                lambda row: f"conn_source_{row['index']}_{box_id}"
                if type(row.id) is float
                else row.id,
                axis=1,
            )

            # Connect targets
            network = snkit.network.link_nodes_to_nearest_edge(
                network,
                lambda node, edge: (
                    (node.asset_type == "target")
                    and (node_to_edge_distance(node.geometry, edge.geometry, geod) < edge_limit)
                )
            )
            network.nodes.loc[network.nodes.id.isnull(), "asset_type"] = "conn_target"

            network.nodes["id"] = network.nodes.reset_index().apply(
                lambda row: f"conn_target_{row['index']}_{box_id}"
                if type(row.id) is float
                else row.id,
                axis=1,
            )

    # Add nodes at line endpoints
    # including where edges have been clipped to slice bbox
    network = snkit.network.add_endpoints(network)

    network.nodes.loc[network.nodes.id.isnull(), "asset_type"] = "intermediate"

    network.nodes["id"] = network.nodes.reset_index().apply(
        lambda row: f"intermediate_{row['index']}_{box_id}"
        if type(row.id) is float
        else row.id,
        axis=1,
    )

    # join econometric and demographic data to network nodes dataframe
    network.nodes = pd.merge(network.nodes, targets[["id", "gdp", "population"]], on="id", how="outer")

    # overwrite ids to int type (grid_disruption.py needs fast pandas isin ops)
    network.nodes["id"] = range(len(network.nodes))
    network.edges["id"] = range(len(network.edges))

    # Add from/to ids
    network = snkit.network.add_topology(network, id_col="id")

    network.edges["box_id"] = box_id
    network.nodes["box_id"] = box_id

    if "component_id" not in network.nodes.columns:
        network = snkit.network.add_component_ids(network)

    network.nodes = allocate_power_to_targets(network.nodes, "gdp")

    # check power has been allocated appropriately
    assert network.nodes["power_mw"].sum() < 1E-6

    network.edges.to_parquet(edges_path)
    network.nodes.to_parquet(nodes_path)
