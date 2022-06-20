#!/usr/bin/env python
# coding: utf-8
"""Read OSM geoparquet, create network, clean it, write out as geopackage."""

import logging
import sys
from typing import Any, Callable
import warnings

import geopandas as gpd

import snkit


def strip_prefix(s: str, prefix: str = 'tag_') -> str:
    """Remove a string prefix if the prefix is present."""

    if s.startswith(prefix):
        return s[len(prefix):]
    else:
        return s


def strip_suffix(s: str, suffix: str = '_link') -> str:
    """Remove a string suffix if the suffix is present."""

    if s.endswith(suffix):
        return s[:len(s) - len(suffix)]
    else:
        return s


def cast(x: Any, *, casting_function: Callable, nullable: bool) -> Any:
    """
    Attempt to recast value with provided function. If not possible and
    nullable is true, return None. Else, raise casting error.

    N.B. Empty string is not considered a nullable value.

    Args:
        x (Any): Value to cast
        casting_function (Callable): Function to cast with
        nullable (bool): Whether cast value can be None

    Returns:
        Any: Recast value
    """

    try:
        new_value = casting_function(x)

        if new_value is None and not nullable:
            raise TypeError(f"{new_value=} with {casting_function=} and {nullable=}")
        else:
            return new_value

    except (ValueError, TypeError) as casting_error:
        # ValueError in case of e.g. x="50 mph"
        # TypeError in case of e.g. x=None

        if nullable:
            return None
        else:
            raise ValueError("Couldn't recast to non-nullable value") from casting_error


def clean_OSM_ways(ways: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Check and clean OpenStreetMap input data

    Args:
        ways (gpd.GeoDataFrame): Table of ways created by osmium and osm_to_pq.py

    Returns:
        gpd.GeoDataFrame: Cleaned table
    """

    # check we have exactly the columns we expect
    necessary_columns = {
        'geometry',
        'tag_highway',
        'tag_surface',
        'tag_bridge',
        'tag_maxspeed',
        'tag_lanes',
    }
    if not necessary_columns.issubset(set(ways.columns.values)):
        raise ValueError(f"{ways.columns.values} does not contain {set(ways.columns.values) - necessary_columns}")

    # check we have only linestrings
    assert set(ways.geometry.type.values) == {'LineString'}

    # recast data where appropriate
    # be careful applying this to string fields, your no-data values may change
    type_conversion_data = (
        ('tag_maxspeed', float, True),
        ('tag_lanes', int, True),
    )
    for column_name, dtype, nullable in type_conversion_data:
        ways[column_name] = ways[column_name].apply(cast, casting_function=dtype, nullable=nullable)

    # make the bridge tag boolean (on read the values are almost entirely None or 'yes')
    ways.tag_bridge = ways.tag_bridge.apply(bool)

    # turn the <highway_type>_link entries into <highway_type>
    ways.tag_highway = ways.tag_highway.apply(strip_suffix)

    # drop the 'tag_' prefix added during osm.pbf processing
    # first check there won't be a collision with new column names
    ways_renamed = ways.rename(strip_prefix, axis='columns')
    n_cols_renamed: int = len(set(ways_renamed.columns.values))
    if n_cols_renamed == len(set(ways.columns.values)):
        ways = ways_renamed
    else:
        raise ValueError(f"Removing prefix would result in collision: {ways.columns=}")

    return ways


def create_network(
    edges: gpd.GeoDataFrame,
    nodes: gpd.GeoDataFrame = None,
    node_edge_prefix: str = ''
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

    if nodes is not None:
        logging.info('Snapping nodes to edges')
        network = snkit.network.snap_nodes(network)

        logging.info('Dropping duplicate geometries')
        network.nodes = snkit.network.drop_duplicate_geometries(network.nodes)

        logging.info('Splitting edges at nodes')
        network = snkit.network.split_edges_at_nodes(network)

    logging.info('Adding endpoints')
    network = snkit.network.add_endpoints(network)

    logging.info('Splitting edges at nodes')
    network = snkit.network.split_edges_at_nodes(network)

    logging.info('Renaming nodes and edges')
    network = snkit.network.add_ids(
        network,
        edge_prefix=f"{node_edge_prefix}e",
        node_prefix=f"{node_edge_prefix}n"
    )

    logging.info('Creating network topology')
    network = snkit.network.add_topology(network, id_col='id')

    network.edges.rename(
        columns={
            'from_id': 'from_node_id',
            'to_id': 'to_node_id',
            'id': 'edge_id'
        },
        inplace=True
    )
    network.nodes.rename(columns={'id': 'node_id'}, inplace=True)

    return network


if __name__ == '__main__':
    try:
        osm_geoparquet_path = snakemake.input[0]
        nodes_output_path = snakemake.output[0]
        edges_output_path = snakemake.output[1]
        network_type = snakemake.config["transport_type"]  # used to label edge IDs and index transport tariffs
        osm_epsg = snakemake.config["osm_epsg"]
    except NameError:
        # If "snakemake" doesn't exist then must be running from the
        # command line.
        osm_geoparquet_path, nodes_output_path, edges_output_path, transport_type, osm_epsg = sys.argv[1:]
        # osm_geoparquet_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0.geoparquet
        # nodes_output_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_nodes.geoparquet
        # edges_output_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_edges.geoparquet
        # transport_type = road
        # osm_epsg = 4326

    osm_epsg = int(osm_epsg)

    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')

    network = create_network(
        edges=clean_OSM_ways(gpd.read_parquet(osm_geoparquet_path)),
        nodes=None,
        node_edge_prefix=network_type
    )
    # manually set crs using geopandas rather than snkit to avoid 'init' style proj crs
    # and permit successful CRS deserializiation and methods such as edges.crs.to_epsg()
    network.edges.set_crs(epsg=osm_epsg, inplace=True)
    network.nodes.set_crs(epsg=osm_epsg, inplace=True)

    logging.info('Writing network to disk')
    network.edges.to_parquet(edges_output_path)
    network.nodes.to_parquet(nodes_output_path)

    logging.info('Done creating network')
