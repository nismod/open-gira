#!/usr/bin/env python
# coding: utf-8
"""Read OSM geoparquet, create network, clean it, write out as geopackage."""

import logging
import sys
from typing import Tuple
import warnings

import geopandas as gpd
import pandas as pd

import utils
from create_network import create_network


def get_rehab_costs(row: pd.Series, rehab_costs: pd.DataFrame) -> Tuple[float, float, str]:
    """
    Determine the cost of rehabilitation for a given rail segment (row).

    Args:
        row (pd.Series): Road segment
        rehab_costs: (pd.DataFrame): Table of rehabilitation costs for various rail types

    Returns:
        Tuple[float, float, str]: Minimum cost, maximum cost, units of cost
    """

    # bridge should be a boolean type after data cleaning step
    if row.bridge:
        asset_type = "bridge"
    else:
        asset_type = "rail"

    data: pd.Series = rehab_costs[rehab_costs["asset_type"] == asset_type].squeeze()

    return data.cost_min, data.cost_max, data.cost_unit


if __name__ == "__main__":
    try:
        osm_edges_path = snakemake.input["edges"]
        osm_nodes_path = snakemake.input["nodes"]
        administrative_data_path = snakemake.input["admin"]
        nodes_output_path = snakemake.output["nodes"]
        edges_output_path = snakemake.output["edges"]
        rehabilitation_costs_path = snakemake.config["transport"]["rehabilitation_costs_path"]
        osm_epsg = snakemake.config["osm_epsg"]
    except NameError:
        # If "snakemake" doesn't exist then must be running from the
        # command line.
        (
            osm_edges_path,
            osm_nodes_path,
            administrative_data_path,
            nodes_output_path,
            edges_output_path,
            rehabilitation_costs_path,
            transport_costs_path,
            flow_cost_time_factor,
            osm_epsg,
        ) = sys.argv[1:]
        # osm_edges_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0.geoparquet
        # osm_nodes_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0.geoparquet
        # administrative_data_path = ../../results/input/admin-boundaries/gadm36_levels.gpkg
        # nodes_output_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_nodes.geoparquet
        # edges_output_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_edges.geoparquet
        # rehabilitation_costs_path = ../../bundled_data/rehabilitation.xlsx
        # osm_epsg = 4326

    osm_epsg = int(osm_epsg)

    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    # read edges
    try:
        edges = gpd.read_parquet(osm_edges_path)
    except ValueError as error:
        # if the input parquet file does not contain a geometry column, geopandas
        # will raise a ValueError rather than try to procede
        logging.info(f"{error}\n" "writing empty files and skipping processing...")
        utils.write_empty_frames(edges_output_path, nodes_output_path)
        sys.exit(0)  # exit gracefully so snakemake will continue

    # read nodes
    try:
        nodes = gpd.read_parquet(osm_nodes_path)
    except ValueError as error:
        # if the input parquet file does not contain a geometry column, geopandas
        # will raise a ValueError rather than try to procede
        logging.info(f"{error}\n" "no nodes from OSM to process...")
        nodes = None

    # if present, filter nodes to stations
    if nodes is not None and not nodes.empty:
        nodes = nodes.loc[nodes.tag_railway == 'station', :]

    network = create_network(edges=edges, nodes=nodes)
    logging.info(
        f"Network contains {len(network.edges)} edges and {len(network.nodes)} nodes"
    )

    # boolean station field
    network.nodes['station'] = network.nodes.tag_railway == 'station'

    # boolean bridge field
    network.edges['bridge'] = utils.str_to_bool(network.edges['tag_bridge'])

    # manually set crs using geopandas rather than snkit to avoid 'init' style proj crs
    # and permit successful CRS deserializiation and methods such as edges.crs.to_epsg()
    network.edges.set_crs(epsg=osm_epsg, inplace=True)
    network.nodes.set_crs(epsg=osm_epsg, inplace=True)

    logging.info("Annotating network with administrative data")
    network = utils.annotate_country(
        network,
        utils.get_administrative_data(administrative_data_path, to_epsg=osm_epsg),
        osm_epsg,
    )

    logging.info("Annotating network with rehabilitation costs")
    network = utils.annotate_rehabilitation_costs(
        network,
        pd.read_excel(rehabilitation_costs_path, sheet_name="rail"),
        get_rehab_costs
    )

    # TODO: add tariffs to edges for flow routing

    # TODO: drop superfluous columns (e.g. OSM tags)

    # TODO: add asset_categories for direct damage estimation

    logging.info("Writing network to disk")
    network.edges.to_parquet(edges_output_path)
    network.nodes.to_parquet(nodes_output_path)

    logging.info("Done creating network")
