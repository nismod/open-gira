#!/usr/bin/env python
# coding: utf-8
"""Read network geopackage, add flow and rehabiliation data, write out."""

import logging
import sys
from typing import Tuple
import warnings

import geopandas as gpd
import pandas as pd
from pyproj import Geod

import snkit

from annotate_road_network import (
    get_administrative_data,
    annotate_country,
    str_to_bool,
    drop_tag_prefix,
)


def clean_edges(edges: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Check and clean OpenStreetMap input data

    Args:
        edges (gpd.GeoDataFrame): Table of edges created by osmium and osm_to_pq.py

    Returns:
        gpd.GeoDataFrame: Cleaned table
    """

    # make the bridge tag boolean
    if "tag_bridge" in edges.columns:
        # coerce our None values into empty strings
        edges.loc[edges['tag_bridge'].isnull(), 'tag_bridge'] = ''
        # these values are mostly a type of bridge construction, 'yes', 'no', or empty string
        edges.tag_bridge = edges.tag_bridge.apply(str_to_bool)

    edges = drop_tag_prefix(edges)

    return edges


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


def annotate_rehabilitation_costs(
    network: snkit.network.Network, rehab_costs: pd.DataFrame
) -> snkit.network.Network:

    # lookup costs
    network.edges["rehab_costs"] = network.edges.apply(
        get_rehab_costs, axis=1, args=(rehab_costs,)
    )

    # unpack results into 3 columns
    network.edges[
        ["rehab_cost_min", "rehab_cost_max", "rehab_cost_unit"]
    ] = network.edges["rehab_costs"].apply(pd.Series)
    network.edges.drop(["rehab_costs"], axis=1, inplace=True)

    return network


if __name__ == "__main__":
    try:
        nodes_path = snakemake.input["nodes"]
        edges_path = snakemake.input["edges"]
        administrative_data_path = snakemake.input["admin"]
        output_nodes_path = snakemake.output["nodes"]
        output_edges_path = snakemake.output["edges"]
        rehabilitation_costs_path = snakemake.config["transport"]["rehabilitation_costs_path"]
        transport_costs_path = snakemake.config["transport"]["tariff_costs_path"]
        flow_cost_time_factor = snakemake.config["transport"]["rail"]["flow_cost_time_factor"]
        osm_epsg = snakemake.config["osm_epsg"]

    except NameError:
        # If "snakemake" doesn't exist then must be running from the
        # command line.
        (
            nodes_path,
            edges_path,
            administrative_data_path,
            output_nodes_path,
            output_edges_path,
            road_speeds_path,
            rehabilitation_costs_path,
            transport_costs_path,
            flow_cost_time_factor,
            osm_epsg,
        ) = sys.argv[1:]
        # nodes_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_nodes.geoparquet
        # edges_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_edges.geoparquet
        # administrative_data_path = ../../results/input/admin-boundaries/gadm36_levels.gpkg
        # output_nodes_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_nodes_annotated.geoparquet
        # output_edges_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_edges_annotated.geoparquet
        # road_speeds_path = ../../bundled_data/global_road_speeds.xlsx
        # rehabilitation_costs_path = ../../bundled_data/rehabilitation_costs.xlsx
        # transport_costs_path = ../../bundled_data/transport_costs.csv
        # flow_cost_time_factor = 0.49
        # osm_epsg = 4326

    # cast script arguments to relevant types where necessary
    flow_cost_time_factor = float(flow_cost_time_factor)
    osm_epsg = int(osm_epsg)

    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    logging.info("Reading network from disk")
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")
    try:
        edges = gpd.read_parquet(edges_path)
        nodes = gpd.read_parquet(nodes_path)
    except ValueError as error:
        # if the input parquet file does not contain a geometry column, geopandas
        # will raise a ValueError rather than try to procede
        logging.info(f"{error}\n" "writing empty files and skipping processing...")

        # snakemake requires that output files exist though, so write empty ones
        empty_gdf = gpd.GeoDataFrame([])
        empty_gdf.to_parquet(output_edges_path)
        empty_gdf.to_parquet(output_nodes_path)
        sys.exit(0)  # exit gracefully so snakemake will continue

    annotated_network = snkit.network.Network(
        edges=clean_edges(gpd.read_parquet(edges_path)),
        nodes=gpd.read_parquet(nodes_path),
    )
    logging.info(
        f"Network contains {len(annotated_network.edges)} edges and {len(annotated_network.nodes)} nodes"
    )

    logging.info("Annotating network with administrative data")
    annotated_network = annotate_country(
        annotated_network,
        get_administrative_data(administrative_data_path, to_epsg=osm_epsg),
        osm_epsg,
    )

    logging.info("Annotating network with rehabilitation costs")
    annotated_network = annotate_rehabilitation_costs(
        annotated_network,
        pd.read_excel(rehabilitation_costs_path, sheet_name="rail"),
    )

    # TODO: add tariffs to edges for flow routing

    # TODO: drop superfluous columns (e.g. OSM tags)

    # TODO: add asset_categories for direct damage estimation

    logging.info("Writing network to disk")
    annotated_network.edges.to_parquet(output_edges_path)
    annotated_network.nodes.to_parquet(output_nodes_path)

    logging.info("Done annotating network with flow and rehabilitation data")
