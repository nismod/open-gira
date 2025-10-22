#!/usr/bin/env python
# coding: utf-8
"""
Read OSM geoparquet, create network, clean it, write out as geopackage.
"""

import logging
import sys
import warnings
from typing import Tuple

import geopandas as gpd
import pandas as pd
import snkit

from utils import annotate_country, cast, strip_suffix
from open_gira.admin import get_administrative_data
from open_gira.assets import RoadAssets
from open_gira.io import write_empty_frames
from open_gira.network_creation import create_network
from open_gira.utils import str_to_bool


def clean_edges(edges: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Check and clean OpenStreetMap input data

    Args:
        edges (gpd.GeoDataFrame): Table of edges created by osmium and osm_to_pq.py

    Returns:
        gpd.GeoDataFrame: Cleaned table
    """

    # recast data where appropriate
    # be careful applying this to string fields, your no-data values may change
    type_conversion_data = (
        ("tag_maxspeed", float, True),
        ("tag_lanes", int, True),
    )
    for column_name, dtype, nullable in type_conversion_data:
        if column_name in edges.columns:
            edges[column_name] = edges[column_name].apply(
                cast, casting_function=dtype, nullable=nullable
            )

    if "tag_highway" in edges.columns:
        # None -> empty string
        edges.loc[edges["tag_highway"].isnull(), "tag_highway"] = ""
        # turn the <highway_type>_link entries into <highway_type>
        edges.tag_highway = edges.tag_highway.apply(strip_suffix)

    # boolean bridge field from tag_bridges
    if "tag_bridge" in edges.columns:
        edges["bridge"] = str_to_bool(edges["tag_bridge"])

    return edges


def get_road_surface(row: pd.Series) -> str:
    """Clean and infer road surface category

    Given a series with 'surface' and 'highway' labels:
    - infer missing surface categories from highway class
    - reclassify paved/unpaved to asphalt/gravel
    - pass any existing surface category unmodified

    N.B. There are several surface categories not considered in this function.
    Here are the major roads recorded for OSM in Tanzania as of June 2022:

    (Pdb) df.tag_surface.value_counts()
    unpaved           2521
    paved             2033
    asphalt           1355
    ground             108
    gravel              38
    compacted           23
    dirt                19
    concrete             5
    concrete:lanes       4
    sand                 2
    fine_gravel          1

    Args:
        row: Must have surface (nullable) and highway attributes.

    Returns:
        surface category string
    """
    if not row.tag_surface:
        if row.tag_highway in {"motorway", "trunk", "primary"}:
            return "asphalt"
        else:
            return None
    elif row.tag_surface == "paved":
        return "asphalt"
    elif row.tag_surface == "unpaved":
        return "gravel"
    else:
        return row.tag_surface


def get_road_is_paved(surface: pd.Series, surface_to_paved: pd.DataFrame) -> pd.Series:
    """Given a series of "surface" categories as strings, infer paved status
    (boolean), default False"""
    dict_ = surface_to_paved.set_index("value").to_dict(orient="dict")["paved"]
    return surface.map(dict_).fillna(False)


def get_road_lanes(row: pd.Series) -> int:
    """
    Given a series characterising a road segment, return a value for the number
    of lanes.

    Args:
        row (pd.Series): Must have have a `lanes` and `highway` labels.

    Returns:
        int: Number of lanes
    """

    try:
        lanes = int(row.tag_lanes)
        if lanes < 1:
            return 1
        else:
            return lanes
    except (ValueError, TypeError):
        # couldn't cast row.lanes into an integer
        # instead guess at lanes from highway classification
        if row.tag_highway in ("motorway", "trunk", "primary"):
            return 2
        else:
            return 1


def annotate_condition(
    network: snkit.network.Network, highway_surface_mapping: pd.DataFrame
) -> snkit.network.Network:
    # infer material type from 'surface' and 'highway' columns
    network.edges["material"] = network.edges.apply(get_road_surface, axis=1)

    # categories materials as paved (true/false)
    network.edges["paved"] = get_road_is_paved(
        network.edges.material, highway_surface_mapping
    )

    # add number of lanes
    network.edges["lanes"] = network.edges.apply(lambda x: get_road_lanes(x), axis=1)

    return network


if __name__ == "__main__":

    osm_edges_path = snakemake.input["edges"]
    osm_nodes_path = snakemake.input["nodes"]
    administrative_data_path = snakemake.input["admin"]
    highway_surface_mapping_path = snakemake.input["highway_surface_mapping"]
    dataset_name = snakemake.wildcards.DATASET
    nodes_output_path = snakemake.output["nodes"]
    edges_output_path = snakemake.output["edges"]
    slice_number = int(snakemake.params["slice_number"])

    osm_epsg = 4326

    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    edges = gpd.read_parquet(osm_edges_path)

    if edges.empty is True:
        write_empty_frames(edges_output_path, nodes_output_path)
        sys.exit(0)  # exit gracefully so snakemake will continue

    # osm_to_pq.py creates these columns but we're not using them, so discard
    edges = edges.drop(
        [
            col
            for col in edges.columns
            if col.startswith("start_node_") or col.startswith("end_node_")
        ],
        axis="columns",
    )

    # for roads we do not currently use any nodes extracted from OSM (osm_nodes_path)
    logging.info("Creating road network")
    network = create_network(
        edges=clean_edges(edges), nodes=None, id_prefix=f"{dataset_name}_{slice_number}"
    )
    logging.info(
        f"Network contains {len(network.edges)} edges and {len(network.nodes)} nodes"
    )

    # manually set crs using geopandas rather than snkit to avoid 'init' style proj crs
    # and permit successful CRS deserializiation and methods such as edges.crs.to_epsg()
    network.edges.set_crs(epsg=osm_epsg, inplace=True)
    network.nodes.set_crs(epsg=osm_epsg, inplace=True)

    logging.info("Annotating network with administrative data")
    network = annotate_country(
        network,
        get_administrative_data(administrative_data_path, to_epsg=osm_epsg),
    )

    logging.info("Annotating network with road type and condition data")
    highway_surface_mapping = pd.read_csv(
        highway_surface_mapping_path, usecols=["value", "paved"], comment="#"
    )
    network = annotate_condition(network, highway_surface_mapping)

    # select and label assets with their type
    # the asset_type is used to later select a damage curve
    # note that order is important here, if an edge is paved, motorway and a bridge, it will be tagged as a bridge only
    network.edges.loc[network.edges.paved == False, "asset_type"] = RoadAssets.UNPAVED
    network.edges.loc[network.edges.paved == True, "asset_type"] = RoadAssets.PAVED
    network.edges.loc[network.edges.tag_highway == "unclassified", "asset_type"] = (
        RoadAssets.UNCLASSIFIED
    )
    network.edges.loc[network.edges.tag_highway == "residential", "asset_type"] = (
        RoadAssets.RESIDENTIAL
    )
    network.edges.loc[network.edges.tag_highway == "tertiary", "asset_type"] = (
        RoadAssets.TERTIARY
    )
    network.edges.loc[network.edges.tag_highway == "secondary", "asset_type"] = (
        RoadAssets.SECONDARY
    )
    network.edges.loc[network.edges.tag_highway == "primary", "asset_type"] = (
        RoadAssets.PRIMARY
    )
    network.edges.loc[network.edges.tag_highway == "trunk", "asset_type"] = (
        RoadAssets.TRUNK
    )
    network.edges.loc[network.edges.tag_highway == "motorway", "asset_type"] = (
        RoadAssets.MOTORWAY
    )
    network.edges.loc[network.edges.bridge == True, "asset_type"] = RoadAssets.BRIDGE

    logging.info("Writing network to disk")
    network.edges.to_parquet(edges_output_path)
    network.nodes.to_parquet(nodes_output_path)

    logging.info("Done creating network")
