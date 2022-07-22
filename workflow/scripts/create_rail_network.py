#!/usr/bin/env python
# coding: utf-8
"""Read OSM geoparquet, create network, clean it, write out as geopackage."""

import logging
import sys
from typing import Any, Callable
import warnings

import geopandas as gpd

import snkit

from create_road_network import create_network


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

    try:
        nodes = gpd.read_parquet(osm_nodes_path)
    except ValueError as error:
        # if the input parquet file does not contain a geometry column, geopandas
        # will raise a ValueError rather than try to procede
        logging.info(f"{error}\n" "no nodes from OSM to process...")
        nodes = None

    if nodes is not None:
        # we only want the station nodes
        nodes = nodes.loc[nodes.tag_railway == 'station', :]

    network = create_network(edges=edges, nodes=nodes)

    # manually set crs using geopandas rather than snkit to avoid 'init' style proj crs
    # and permit successful CRS deserializiation and methods such as edges.crs.to_epsg()
    network.edges.set_crs(epsg=osm_epsg, inplace=True)
    network.nodes.set_crs(epsg=osm_epsg, inplace=True)

    logging.info("Writing network to disk")
    network.edges.to_parquet(edges_output_path)
    network.nodes.to_parquet(nodes_output_path)

    logging.info("Done creating network")
