"""
Creates the network from the plants and targets data
"""
import os

import geopandas
import pandas
import snkit
import snkit.network
from pyproj import Geod


if __name__ == "__main__":
    output_dir = snakemake.config["output_dir"]  # type: ignore
    box_id = snakemake.wildcards.BOX  # type: ignore
    plants_path = snakemake.input.plants  # type: ignore
    targets_path = snakemake.input.targets  # type: ignore
    gridfinder_path = snakemake.input.gridfinder  # type: ignore
    nodes_path = snakemake.output.nodes  # type: ignore
    edges_path = snakemake.output.edges  # type: ignore

    # Nodes
    plants = geopandas.read_parquet(plants_path)
    targets = geopandas.read_parquet(targets_path)

    plants["id"] = [f"source_{i}_{box_id}" for i in range(len(plants))]
    plants = plants[["id", "source_id", "capacity_mw", "type", "geometry"]]

    targets["id"] = [f"target_{i}_{box_id}" for i in range(len(targets))]
    target_cols = list(targets.columns)

    target_nodes = targets[["id", "type", "geometry"]].copy()
    target_nodes.geometry = target_nodes.geometry.centroid

    nodes = geopandas.GeoDataFrame(
        pandas.concat([plants, target_nodes], ignore_index=True),
        crs=plants.crs)

    # Edges
    edges = geopandas.read_parquet(gridfinder_path)

    edges["type"] = "transmission"
    edges["id"] = edges.reset_index()["index"].apply(
        lambda i: f"edge_{i}_{box_id}"
    )
    edges = edges[["id", "source_id", "type", "geometry", "source"]]

    # Process network
    network = snkit.network.Network(nodes, edges)

    if len(edges) > 0:
        network = snkit.network.split_multilinestrings(network)

        geod = Geod(ellps="WGS84")
        edge_limit = 20_000  # meters - TODO test limit

        # Connect power plants
        network = snkit.network.link_nodes_to_nearest_edge(
            network,
            lambda node, edge: node.type == "source"
            and geod.geometry_length(edge.geometry) < edge_limit,
        )
        network.nodes.loc[network.nodes.id.isnull(), "type"] = "conn_source"

        network.nodes["id"] = network.nodes.reset_index().apply(
            lambda row: f"conn_source_{row['index']}_{box_id}"
            if type(row.id) is float
            else row.id,
            axis=1,
        )

        # Connect targets
        network = snkit.network.link_nodes_to_nearest_edge(
            network,
            lambda node, edge: node.type == "target"
            and geod.geometry_length(edge.geometry) < edge_limit,
        )
        network.nodes.loc[network.nodes.id.isnull(), "type"] = "conn_target"

        network.nodes["id"] = network.nodes.reset_index().apply(
            lambda row: f"conn_target_{row['index']}_{box_id}"
            if type(row.id) is float
            else row.id,
            axis=1,
        )

    # Add nodes at line endpoints
    network = snkit.network.add_endpoints(network)

    network.nodes.loc[network.nodes.id.isnull(), "type"] = "intermediate"

    network.nodes["id"] = network.nodes.reset_index().apply(
        lambda row: f"intermediate_{row['index']}_{box_id}"
        if type(row.id) is float
        else row.id,
        axis=1,
    )

    # Add from/to ids
    network = snkit.network.add_topology(network, id_col="id")

    network.edges["box_id"] = box_id
    network.nodes["box_id"] = box_id

    network.edges.to_parquet(edges_path)
    network.nodes.to_parquet(nodes_path)
