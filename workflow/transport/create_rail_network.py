#!/usr/bin/env python
# coding: utf-8
"""
Read OSM geoparquet, create network, clean it, write out as geopackage.
"""

import logging
import sys
import warnings

import geopandas as gpd
import pandas as pd
import snkit

from utils import annotate_country, get_administrative_data
from open_gira.assets import RailAssets
from open_gira.io import write_empty_frames
from open_gira.network import create_network
from open_gira.utils import str_to_bool


def annotate_rehabilitation_costs(
    network: snkit.network.Network, rehab_costs: pd.DataFrame
) -> snkit.network.Network:

    def get_rehab_costs(row: pd.Series, rehab_costs: pd.DataFrame) -> float:
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

        return float(rehab_costs.loc[rehab_costs["asset_type"] == asset_type, "rehab_cost_USD_per_km"])

    # lookup costs
    network.edges["rehab_cost_USD_per_km"] = network.edges.apply(
        get_rehab_costs, axis=1, args=(rehab_costs,)
    )

    return network


if __name__ == "__main__":

    try:
        osm_edges_path = snakemake.input["edges"]  # type: ignore
        osm_nodes_path = snakemake.input["nodes"]  # type: ignore
        administrative_data_path = snakemake.input["admin"]  # type: ignore
        nodes_output_path = snakemake.output["nodes"]  # type: ignore
        edges_output_path = snakemake.output["edges"]  # type: ignore
        slice_number = snakemake.params["slice_number"]  # type: ignore
        rehabilitation_costs_path = snakemake.config["transport"]["rehabilitation_costs_path"]  # type: ignore

    except NameError:
        # If "snakemake" doesn't exist then must be running from the
        # command line.
        (
            osm_edges_path,
            osm_nodes_path,
            administrative_data_path,
            nodes_output_path,
            edges_output_path,
            slice_number,
            rehabilitation_costs_path,
            transport_costs_path,
            flow_cost_time_factor,
        ) = sys.argv[1:]

        # osm_edges_path = ../../results/geoparquet/tanzania-latest_filter-rail/slice-0_edges.geoparquet
        # osm_nodes_path = ../../results/geoparquet/tanzania-latest_filter-rail/slice-0_nodes.geoparquet
        # administrative_data_path = ../../results/input/admin-boundaries/gadm36_levels.gpkg
        # nodes_output_path = ../../results/geoparquet/tanzania-latest_filter-rail/slice-0_nodes.geoparquet
        # edges_output_path = ../../results/geoparquet/tanzania-latest_filter-rail/slice-0_edges.geoparquet
        # slice_number = 0
        # rehabilitation_costs_path = ../../bundled_data/rehabilitation.xlsx

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

    # select and label assets with their type
    network.nodes.loc[network.nodes.tag_railway == 'station', 'asset_type'] = RailAssets.STATION
    network.edges.loc[str_to_bool(network.edges['tag_bridge']), 'asset_type'] = RailAssets.BRIDGE
    network.edges.loc[network.edges.tag_railway == 'rail', 'asset_type'] = RailAssets.RAILWAY

    # boolean station field
    network.nodes['station'] = network.nodes.tag_railway == 'station'

    # boolean bridge field
    network.edges['bridge'] = str_to_bool(network.edges['tag_bridge'])

    # manually set crs using geopandas rather than snkit to avoid 'init' style proj crs
    # and permit successful CRS deserializiation and methods such as edges.crs.to_epsg()
    network.edges.set_crs(epsg=osm_epsg, inplace=True)
    network.nodes.set_crs(epsg=osm_epsg, inplace=True)

    logging.info("Annotating network with administrative data")
    network = annotate_country(
        network,
        get_administrative_data(administrative_data_path, to_epsg=osm_epsg),
    )

    logging.info("Annotate network with rehabilitation costs")
    network = annotate_rehabilitation_costs(
        network,
        pd.read_excel(rehabilitation_costs_path, sheet_name="rail")
    )

    logging.info("Writing network to disk")
    network.edges.to_parquet(edges_output_path)
    network.nodes.to_parquet(nodes_output_path)

    logging.info("Done creating network")
