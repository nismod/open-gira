#!/usr/bin/env python
# coding: utf-8
"""Read OSM geoparquet, create network, clean it, write out as geopackage."""

import logging
import sys
from typing import Any, Callable
import warnings

import geopandas as gpd

import snkit


def create_network(
    edges: gpd.GeoDataFrame, nodes: gpd.GeoDataFrame = None, node_edge_prefix: str = ""
) -> snkit.network.Network:
    """
    Create snkit network from edges and (optional) nodes and clean the result.

    Arguments:
        edges (gpd.GeoDataFrame): Expected to contain geometry column of linestrings
        nodes (gpd.GeoDataFrame): Optional nodes to include. snkit will try to snap to edges
        node_edge_prefix (str): Returned network will have nodes and edge IDs prefixed by this string

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
    network = snkit.network.add_ids(network)

    logging.info("Creating network topology")
    network = snkit.network.add_topology(network, id_col="id")

    network.edges.rename(
        columns={"from_id": "from_node_id", "to_id": "to_node_id", "id": "edge_id"},
        inplace=True,
    )
    network.nodes.rename(columns={"id": "node_id"}, inplace=True)

    return network


if __name__ == "__main__":
    try:
        osm_edges_path = snakemake.input["edges"]
        osm_nodes_path = snakemake.input["nodes"]
        nodes_output_path = snakemake.output["nodes"]
        edges_output_path = snakemake.output["edges"]
        osm_epsg = snakemake.config["osm_epsg"]
    except NameError:
        # If "snakemake" doesn't exist then must be running from the
        # command line.
        (
            osm_edges_path,
            osm_nodes_path,
            nodes_output_path,
            edges_output_path,
            osm_epsg,
        ) = sys.argv[1:]
        # osm_edges_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0.geoparquet
        # osm_nodes_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0.geoparquet
        # nodes_output_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_nodes.geoparquet
        # edges_output_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_edges.geoparquet
        # osm_epsg = 4326

    osm_epsg = int(osm_epsg)

    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    try:
        edges = gpd.read_parquet(osm_edges_path)
    except ValueError as error:
        # if the input parquet file does not contain a geometry column, geopandas
        # will raise a ValueError rather than try to procede
        logging.info(f"{error}\n" "writing empty files and skipping processing...")

        # snakemake requires that output files exist though, so write empty ones
        empty_gdf = gpd.GeoDataFrame([])
        empty_gdf.to_parquet(edges_output_path)
        empty_gdf.to_parquet(nodes_output_path)
        sys.exit(0)  # exit gracefully so snakemake will continue

    # for roads we do not currently use any nodes extracted from OSM (osm_nodes_path)
    network = create_network(edges=edges, nodes=None)

    # manually set crs using geopandas rather than snkit to avoid 'init' style proj crs
    # and permit successful CRS deserializiation and methods such as edges.crs.to_epsg()
    network.edges.set_crs(epsg=osm_epsg, inplace=True)
    network.nodes.set_crs(epsg=osm_epsg, inplace=True)

    logging.info("Writing network to disk")
    network.edges.to_parquet(edges_output_path)
    network.nodes.to_parquet(nodes_output_path)

    logging.info("Done creating network")
