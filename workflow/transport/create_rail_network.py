#!/usr/bin/env python
# coding: utf-8
"""
Read OSM geoparquet, create network, clean it, write out as geopackage.
"""

import logging
import sys
import warnings

import geopandas as gpd

from utils import annotate_country, get_administrative_data
from open_gira.assets import RailAssets
from open_gira.io import write_empty_frames
from open_gira.network import create_network
from open_gira.utils import str_to_bool


if __name__ == "__main__":

    osm_edges_path = snakemake.input["edges"]  # type: ignore
    osm_nodes_path = snakemake.input["nodes"]  # type: ignore
    administrative_data_path = snakemake.input["admin"]  # type: ignore
    nodes_output_path = snakemake.output["nodes"]  # type: ignore
    edges_output_path = snakemake.output["edges"]  # type: ignore
    slice_number = int(snakemake.params["slice_number"])  # type: ignore

    slice_number = int(slice_number)
    osm_epsg = 4326

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    # read edges
    edges = gpd.read_parquet(osm_edges_path)
    if edges.empty is True:
        write_empty_frames(edges_output_path, nodes_output_path)
        sys.exit(0)  # exit gracefully so snakemake will continue

    # read nodes
    nodes = gpd.read_parquet(osm_nodes_path)
    if nodes.empty is True:
        nodes = None

    # osm_to_pq.py creates these columns but we're not using them, so discard
    edges = edges.drop(
        [col for col in edges.columns if col.startswith("start_node_") or col.startswith("end_node_")],
        axis="columns"
    )

    # if present, filter nodes to stations
    if nodes is not None and not nodes.empty:
        nodes = nodes.loc[nodes.tag_railway == 'station', :]

    # pass an id_prefix containing the slice number to ensure edges and nodes
    # are uniquely identified across all slices in the network
    network = create_network(edges=edges, nodes=nodes, id_prefix=f"{slice_number}")
    logging.info(
        f"Network contains {len(network.edges)} edges and {len(network.nodes)} nodes"
    )

    # boolean bridge field
    network.edges['bridge'] = str_to_bool(network.edges['tag_bridge'])

    # boolean station field
    network.nodes['station'] = network.nodes.tag_railway == 'station'

    # select and label assets with their type
    # we will use the `asset_type` field to select damage curves
    # bridge overrides railway as asset class, tag last
    network.nodes.loc[network.nodes.station == True, 'asset_type'] = RailAssets.STATION
    network.edges.loc[network.edges.tag_railway == 'rail', 'asset_type'] = RailAssets.RAILWAY
    network.edges.loc[network.edges.bridge == True, 'asset_type'] = RailAssets.BRIDGE

    # manually set crs using geopandas rather than snkit to avoid 'init' style proj crs
    # and permit successful CRS deserializiation and methods such as edges.crs.to_epsg()
    network.edges.set_crs(epsg=osm_epsg, inplace=True)
    network.nodes.set_crs(epsg=osm_epsg, inplace=True)

    logging.info("Annotating network with administrative data")
    network = annotate_country(
        network,
        get_administrative_data(administrative_data_path, to_epsg=osm_epsg),
    )

    logging.info("Writing network to disk")
    network.edges.to_parquet(edges_output_path)
    network.nodes.to_parquet(nodes_output_path)

    logging.info("Done creating network")