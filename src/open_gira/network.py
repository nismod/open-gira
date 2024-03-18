"""
Functionality for creating and manipulating networks.
"""

import logging
import warnings

import geopandas as gpd
import snkit


def create_network(
    edges: gpd.GeoDataFrame,
    nodes: gpd.GeoDataFrame = None,
    id_prefix: str = ""
) -> snkit.network.Network:
    """
    Create snkit network from edges and (optional) nodes and clean the result.

    Arguments:
        edges (gpd.GeoDataFrame): Expected to contain geometry column of linestrings
        nodes (gpd.GeoDataFrame): Optional nodes to include. snkit will try to snap to edges

    Returns:
        snkit.network.Network: Built network
    """

    logging.info("Starting network creation")

    # drop edges with no geometry
    empty_idx = edges.geometry.apply(lambda e: e is None or e.is_empty)
    if empty_idx.sum():
        empty_edges = edges[empty_idx]
        logging.info(f"Found {len(empty_edges)} empty edges.")
        logging.info(empty_edges)
        edges = edges[~empty_idx].copy()

    logging.info("Creating network")
    network = snkit.Network(nodes, edges)

    logging.info("Splitting multilines")
    network = snkit.network.split_multilinestrings(network)

    if nodes is not None and not nodes.empty:
        # check we have only point nodes
        assert set(network.nodes.geometry.type.values) == {"Point"}

        logging.info("Dropping duplicate geometries")
        # silence shapely.ops.split ShapelyDeprecationWarning regarding:
        # shapley.ops.split failure on split by empty geometry collection
        # this is currently caught by snkit.network.split_edge_at_points,
        # but won't be for shapely==2.0
        warnings.filterwarnings(
            "ignore",
            message=(
                ".*GeometryTypeError will derive from ShapelyError "
                "and not TypeError or ValueError in Shapely 2.0*"
            )
        )
        network.nodes = snkit.network.drop_duplicate_geometries(network.nodes)

        logging.info("Snapping nodes to edges")
        network = snkit.network.snap_nodes(network)

    logging.info("Adding endpoints")
    network = snkit.network.add_endpoints(network)

    logging.info("Splitting edges at nodes")
    network = snkit.network.split_edges_at_nodes(network)

    # check we have only linestrings
    assert set(network.edges.geometry.type.values) == {"LineString"}

    logging.info("Renaming nodes and edges")
    network = snkit.network.add_ids(network, edge_prefix=id_prefix, node_prefix=id_prefix)

    logging.info("Creating network topology")
    network = snkit.network.add_topology(network, id_col="id")

    return network
