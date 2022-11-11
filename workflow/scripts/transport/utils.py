#!/usr/bin/env python
# coding: utf-8
"""
Shared functions for creating, cleaning, manipulating and analysing networks.
"""

import logging
import re
from typing import Any, Callable

import geopandas as gpd
import pandas as pd
import snkit


WEB_MERC_EPSG = 3857  # Web Mercator, a projected CRS


def strip_prefix(s: str, prefix: str = "tag_") -> str:
    """Remove a string prefix if the prefix is present."""
    return re.sub(f"^{prefix}", "", s)


def strip_suffix(s: str, suffix: str = "_link") -> str:
    """Remove a string suffix if the suffix is present."""
    return re.sub(f"{suffix}$", "", s)


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


def get_administrative_data(file_path: str, to_epsg: int = None) -> gpd.GeoDataFrame:
    """
    Read administrative data (country ISO, country geometry) from disk

    Arguments:
        file_path (str): Location of file with country data:
            containing an ISO three letter code as 'GID_0' and a geometry as
            'geometry'
        to_epsg (int): EPSG code to project data to

    Returns:
        gpd.GeoDataFrame: Table of country and geometry data with:
            'iso_a3' and 'geometry' columns
    """

    # read file
    gdf = gpd.read_file(file_path)

    # check schema is as expected
    expected_columns = {"GID_0", "geometry"}
    assert expected_columns.issubset(set(gdf.columns.values))

    # reproject if desired
    if to_epsg is not None:
        gdf = gdf.to_crs(epsg=to_epsg)

    # rename these columns first so we don't have to do this twice (to nodes and edges) later
    gdf.rename(columns={"GID_0": "iso_a3"}, inplace=True)

    # subset, sort by iso_a3 and return
    return gdf[["iso_a3", "geometry"]].sort_values(by=["iso_a3"], ascending=True)


def annotate_country(network: snkit.network.Network, countries: gpd.GeoDataFrame) -> snkit.network.Network:
    """
    Label network edges and nodes with their country ISO code

    Arguments:
        network (snkit.network.Network): Network to label with geographic information
            network.edges must have 'from_node_id', 'to_node_id'
            network.nodes must have 'node_id', 'geometry'
        countries (gpd.GeoDataFrame): Table required to contain the following columns:
            'iso_a3', 'geometry'

    Returns
        snkit.network.Network: Labelled network with edges.from_iso_a3, edges.to_iso_a3
            and nodes.iso_a3
    """
    # unpack network for brevity in the following lines
    edges = network.edges
    nodes = network.nodes

    # check we have required inputs
    assert set(edges.columns.values).issuperset({"from_node_id", "to_node_id"})
    assert set(nodes.columns.values).issuperset({"node_id", "geometry"})
    assert set(countries.columns.values).issuperset({"iso_a3", "geometry"})

    # check our inputs have a registered CRS
    assert nodes.crs
    # we don't do any geometry operations on edges
    assert countries.crs

    input_node_crs = nodes.crs

    # shorthand for selecting certain column sets
    starting_node_columns = list(nodes.columns.values)
    desired_node_columns = starting_node_columns + ['iso_a3']

    # spatial join nodes geometries to their containing country, retain only node geometries
    nodes_with_iso_a3 = nodes.sjoin(countries, how="left", predicate="within")
    # drop cruft from merge (i.e. "index_right")
    nodes_with_iso_a3 = nodes_with_iso_a3[desired_node_columns]
    interior_nodes = nodes_with_iso_a3[~nodes_with_iso_a3["iso_a3"].isna()]
    logging.info(f"Found {len(interior_nodes)} nodes that are within a country geometry")

    # for any nodes where sjoin didn't work, drop the iso_a3 and try again with sjoin_nearest
    # reproject to web mercator CRS for this operation
    exterior_nodes = nodes_with_iso_a3[nodes_with_iso_a3["iso_a3"].isna()].drop(["iso_a3"], axis="columns")
    exterior_nodes = exterior_nodes.to_crs(WEB_MERC_EPSG).sjoin_nearest(countries.to_crs(WEB_MERC_EPSG))
    exterior_nodes = exterior_nodes.to_crs(input_node_crs)
    if not exterior_nodes.empty:
        logging.info(
            f"Found {len(exterior_nodes)} nodes external to all country geometries, labelled "
            "these with their nearest country"
        )

    # join the outputs of the two joining processes together
    nodes = pd.concat([interior_nodes, exterior_nodes], axis=0, ignore_index=True)
    nodes = gpd.GeoDataFrame(
        # use the columns we were passed, plus the country code of each node
        nodes[desired_node_columns],
        geometry="geometry",
        crs=input_node_crs
    )

    # set edge.from_node_id from node.node_id and use iso_a3 of from node as edge start
    edges = pd.merge(
        edges,
        nodes[["node_id", "iso_a3"]],
        how="left",
        left_on=["from_node_id"],
        right_on=["node_id"],
    )
    edges.rename(columns={"iso_a3": "from_iso_a3"}, inplace=True)
    edges.drop("node_id", axis=1, inplace=True)

    # set edge.to_node_id from node.node_id and use iso_a3 of from node as edge end
    edges = pd.merge(
        edges,
        nodes[["node_id", "iso_a3"]],
        how="left",
        left_on=["to_node_id"],
        right_on=["node_id"],
    )
    edges.rename(columns={"iso_a3": "to_iso_a3"}, inplace=True)
    edges.drop("node_id", axis=1, inplace=True)

    network.nodes = nodes
    network.edges = edges

    return network
