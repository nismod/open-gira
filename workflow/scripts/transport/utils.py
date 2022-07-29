#!/usr/bin/env python
# coding: utf-8
"""
Shared functions for creating, cleaning, manipulating and analysing networks.
"""


import logging
from typing import Any, Callable

import pandas as pd
import geopandas as gpd

import snkit


WEB_MERC_EPSG = 3857  # Web Mercator, a projected CRS


def write_empty_frames(edges_path: str, nodes_path: str) -> None:
    """
    If we don't have sufficient / good enough input data, write out empty output.

    N.B. Output files must exist for snakemake's sake.
    """
    empty_gdf = gpd.GeoDataFrame([])
    empty_gdf.to_parquet(edges_path)
    empty_gdf.to_parquet(nodes_path)
    return


def strip_prefix(s: str, prefix: str = "tag_") -> str:
    """Remove a string prefix if the prefix is present."""

    if s.startswith(prefix):
        return s[len(prefix):]
    else:
        return s


def strip_suffix(s: str, suffix: str = "_link") -> str:
    """Remove a string suffix if the suffix is present."""
    if not isinstance(s, str):
        raise ValueError

    if s.endswith(suffix):
        return s[: len(s) - len(suffix)]
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


def str_to_bool(series: pd.Series) -> pd.Series:
    """
    Make a stab at turning a series of strings into a boolean series.

    Args:
        series (pd.Series): Input series

    Returns:
        pd.Series: Boolean series of whether or not strings are truthy (by our logic).
    """

    FALSE_VALUES = {'n', 'no', 'false', 'f', ''}

    def str_parser(s: str) -> bool:
        """
        If a string is in the set below, return False, otherwise, return True.
        """

        if not isinstance(s, str):
            raise ValueError(f"{s=} has {type(s)=}, but should be str")

        return s.lower() not in FALSE_VALUES

    # set our null values to the empty string
    new_series = series.copy(deep=True)
    new_series.loc[series.isnull()] = ''

    return new_series.apply(str_parser)


def get_administrative_data(file_path: str, to_epsg: int = None) -> gpd.GeoDataFrame:
    """
    Read administrative data (country ISO, country geometry) from disk

    Arguments:
        file_path (str): Location of file with country data:
            containing an iso_code as 'GID_0' and a geometry as 'geometry'
        to_epsg (int): EPSG code to project data to

    Returns:
        gpd.GeoDataFrame: Table of country and geometry data with:
            'iso_code' and 'geometry' columns
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
    gdf.rename(columns={"GID_0": "iso_code"}, inplace=True)

    # subset, sort by iso_code and return
    return gdf[["iso_code", "geometry"]].sort_values(by=["iso_code"], ascending=True)


def annotate_country(network: snkit.network.Network, countries: gpd.GeoDataFrame) -> snkit.network.Network:
    """
    Label network edges and nodes with their country ISO code

    Arguments:
        network (snkit.network.Network): Network to label with geographic information
            network.edges must have 'from_node_id', 'to_node_id'
            network.nodes must have 'node_id', 'geometry'
        countries (gpd.GeoDataFrame): Table required to contain the following columns:
            'iso_code', 'geometry'

    Returns
        snkit.network.Network: Labelled network with edges.from_iso_code, edges.to_iso_code
            and nodes.iso_code
    """
    # unpack network for brevity in the following lines
    edges = network.edges
    nodes = network.nodes

    # check we have required inputs
    assert set(edges.columns.values).issuperset({"from_node_id", "to_node_id"})
    assert set(nodes.columns.values).issuperset({"node_id", "geometry"})
    assert set(countries.columns.values).issuperset({"iso_code", "geometry"})

    # check our inputs have a registered CRS
    assert nodes.crs
    # we don't do any geometry operations on edges
    assert countries.crs

    input_node_crs = nodes.crs

    # shorthand for selecting certain column sets
    starting_node_columns = list(nodes.columns.values)
    desired_node_columns = starting_node_columns + ['iso_code']

    # spatial join nodes geometries to their containing country, retain only node geometries
    nodes_with_iso_code = nodes.sjoin(countries, how="left", predicate="within")
    # drop cruft from merge (i.e. "index_right")
    nodes_with_iso_code = nodes_with_iso_code[desired_node_columns]
    interior_nodes = nodes_with_iso_code[~nodes_with_iso_code["iso_code"].isna()]
    logging.info(f"Found {len(interior_nodes)} nodes that are within a country geometry")

    # for any nodes where sjoin didn't work, drop the iso_code and try again with sjoin_nearest
    # reproject to web mercator CRS for this operation
    exterior_nodes = nodes_with_iso_code[nodes_with_iso_code["iso_code"].isna()].drop(["iso_code"], axis="columns")
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

    # set edge.from_node_id from node.node_id and use iso_code of from node as edge start
    edges = pd.merge(
        edges,
        nodes[["node_id", "iso_code"]],
        how="left",
        left_on=["from_node_id"],
        right_on=["node_id"],
    )
    edges.rename(columns={"iso_code": "from_iso_code"}, inplace=True)
    edges.drop("node_id", axis=1, inplace=True)

    # set edge.to_node_id from node.node_id and use iso_code of from node as edge end
    edges = pd.merge(
        edges,
        nodes[["node_id", "iso_code"]],
        how="left",
        left_on=["to_node_id"],
        right_on=["node_id"],
    )
    edges.rename(columns={"iso_code": "to_iso_code"}, inplace=True)
    edges.drop("node_id", axis=1, inplace=True)

    network.nodes = nodes
    network.edges = edges

    return network


def annotate_rehabilitation_costs(
    network: snkit.network.Network, rehab_costs: pd.DataFrame, getter: Callable
) -> snkit.network.Network:

    # lookup costs
    network.edges["rehab_costs"] = network.edges.apply(
        getter, axis=1, args=(rehab_costs,)
    )

    # unpack results into 3 columns
    network.edges[
        ["rehab_cost_min", "rehab_cost_max", "rehab_cost_unit"]
    ] = network.edges["rehab_costs"].apply(pd.Series)
    network.edges.drop(["rehab_costs"], axis=1, inplace=True)

    return network
