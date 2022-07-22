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
    annotate_asset_category,
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


def annotate_tariff_flow_costs(
    network: snkit.network.Network,
    transport_tariffs: pd.DataFrame,
    flow_cost_time_factor: float,
) -> snkit.network.Network:
    """
    Add tariff flow costs to network edges.

    Args:
        network (snkit.network.Network): Network to annotate. The network.edges
            must have a 'from_iso' column to merge on.
        transport_tariffs (pd.DataFrame): Table of transport tariffs by country and
            by mode of transport, with 'from_iso3', 'transport', 'cost_km', 'cost_unit'
            and 'cost_scaling' columns
        flow_cost_time_factor (float): A fudge factor that varies (by country?)
            This may well need consuming as location specific data in future.
    Returns:
        snkit.network.Network: Modified network
    """

    # input checking
    expected_columns = {
        "from_iso3",
        "cost_km",
        "cost_unit",
        "cost_scaling",
    }
    if not set(transport_tariffs.columns).issuperset(expected_columns):
        raise ValueError(f"{expected_columns=} for transport_tariffs")

    # rename and subset table
    transport_tariffs.rename(
        columns={"cost_km": "tariff_cost", "cost_unit": "tariff_unit"}, inplace=True
    )

    # merge datasets
    network.edges = pd.merge(
        network.edges,
        transport_tariffs[["from_iso3", "tariff_cost", "tariff_unit"]],
        how="left",
        left_on=["from_iso"],
        right_on=["from_iso3"],
    )
    network.edges["min_tariff"] = network.edges.apply(
        lambda x: float(x.tariff_cost) - (float(x.tariff_cost) * 0.2), axis=1
    )
    network.edges["max_tariff"] = network.edges.apply(
        lambda x: float(x.tariff_cost) + (float(x.tariff_cost) * 0.2), axis=1
    )
    network.edges.drop(["tariff_cost", "from_iso3"], axis=1, inplace=True)

    # assign flow costs
    metres_per_km = 1_000

    # calculate rail segment lengths
    geod = Geod(ellps="WGS84")
    network.edges["length_m"] = network.edges.apply(
        lambda x: float(geod.geometry_length(x.geometry)), axis=1
    )

    # TODO: assign simple transport cost to permit rudimentary flow modelling
    #network.edges["min_flow_cost"] = ?
    #network.edges["max_flow_cost"] = ?
    #network.edges["flow_cost_unit"] = ?

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
        asset_categories = snakemake.config["direct_damage"]["asset_categories"]

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
            asset_categories,
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
        # asset_categories = 'road_paved, road_unpaved, road_bridge'

        # in case of CLI asset_categories argument, remove whitespace and then split on comma
        asset_categories: list[str] = "".join(asset_categories.split()).split(',')

    # cast script arguments to relevant types where necessary
    flow_cost_time_factor = float(flow_cost_time_factor)
    osm_epsg = int(osm_epsg)
    asset_categories = set(asset_categories)

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

    logging.info("Annotating network with tariff and flow costs")
    annotated_network = annotate_tariff_flow_costs(
        annotated_network,
        pd.read_excel(transport_costs_path, sheet_name="rail"),
        flow_cost_time_factor,
    )

    logging.info("Annotating network with asset category")
    annotated_network = annotate_asset_category(annotated_network, asset_categories)

    # TODO: drop superfluous columns (e.g. OSM tags)

    logging.info("Writing network to disk")
    annotated_network.edges.to_parquet(output_edges_path)
    annotated_network.nodes.to_parquet(output_nodes_path)

    logging.info("Done annotating network with flow and rehabilitation data")
