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


def allocate_power(network: snkit.network.Network, target_gdp: gpd.GeoDataFrame) -> snkit.network.Network:
    """
    Allocate total generating capacity to targets (consuming nodes), weighted
    by targets' GDP.

    Given a `network` and the following `target_gdp`:
                        id           gdp
    914    target_914_1030  8.256908e+07
    951    target_951_1030  1.441001e+06
    952    target_952_1030  2.087217e+06
    958    target_958_1030  2.353413e+06
    959    target_959_1030  1.342598e+07
    967    target_967_1030  3.044258e+05

    Here's the returned `network.nodes`:
                        id   source_id  power_mw    asset_type  component_id
    37      source_37_1030  WRI1028006  6.000000        source             6
    956    target_914_1030         NaN -4.848396        target             6
    993    target_951_1030         NaN -0.084615        target             6
    994    target_952_1030         NaN -0.122560        target             6
    1000   target_958_1030         NaN -0.138191        target             6
    1001   target_959_1030         NaN -0.788364        target             6
    1009   target_967_1030         NaN -0.017876        target             6

    Args:
        network: network.nodes object should contain `id`, `power_mw` and
            `asset_type` columns
        targets: should contain `id` and `gdp` columns

    Returns:
        Network with generating capacity allocated to consuming nodes
    """

    if "component_id" not in network.nodes.columns:
        network = snkit.network.add_component_ids(network)

    nodes = pd.merge(network.nodes, target_gdp[["id", "gdp"]], on="id", how="outer")

    # for each component, allocate generation weighted by GDP of targets
    for c_id in set(nodes.component_id):

        c_mask = nodes.component_id == c_id

        c_total_capacity: float = nodes.loc[(nodes.asset_type == "source") & c_mask, "power_mw"].sum()

        c_target_mask = (nodes.asset_type == "target") & c_mask

        c_total_gdp: float = nodes.loc[c_target_mask, "gdp"].sum()

        nodes.loc[c_target_mask, "power_mw"] = \
            -1 * c_total_capacity * (nodes.loc[c_target_mask, "gdp"] / c_total_gdp)

    # check power has been allocated to within a Watt
    assert nodes.power_mw.sum() < 1E-6

    network.nodes = nodes.drop(columns=["gdp"])

    return network


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

    network.edges["id"] = [f"edge_{i}_{box_id}" for i in range(len(network.edges))]

    # Add from/to ids
    network = snkit.network.add_topology(network, id_col="id")

    network.edges["box_id"] = box_id
    network.nodes["box_id"] = box_id

    network = allocate_power(network, targets)

    network.edges.to_parquet(edges_path)
    network.nodes.to_parquet(nodes_path)
