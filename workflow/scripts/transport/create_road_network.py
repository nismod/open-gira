#!/usr/bin/env python
# coding: utf-8
"""Read OSM geoparquet, create network, clean it, write out as geopackage.
"""
import logging
import sys
import warnings
from typing import Tuple

import geopandas as gpd
import pandas as pd
import snkit
from pyproj import Geod

from create_network import create_network
from assets import RoadAssets
from utils import (annotate_country, annotate_rehabilitation_costs, cast,
                    get_administrative_data, str_to_bool, strip_suffix,
                    write_empty_frames)


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
        edges.loc[edges['tag_highway'].isnull(), 'tag_highway'] = ''
        # turn the <highway_type>_link entries into <highway_type>
        edges.tag_highway = edges.tag_highway.apply(strip_suffix)

    # boolean bridge field from tag_bridges
    if "tag_bridge" in edges.columns:
        edges['bridge'] = str_to_bool(edges['tag_bridge'])

    return edges


def get_road_condition(row: pd.Series) -> Tuple[str, str]:
    """
    Given a series with 'surface' and 'highway' labels, infer road:
        - paved status (boolean)
        - surface category from {'asphalt', 'gravel', 'concrete'}

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
        row (pd.Series): Must have surface (nullable) and highway attributes.

    Returns:
        Tuple[str, str]: Categorical condition and surface strings
    """

    if not row.tag_surface:
        if row.tag_highway in {"motorway", "trunk", "primary"}:
            return True, "asphalt"
        else:
            return False, "gravel"
    elif row.tag_surface == "paved":
        return True, "asphalt"
    elif row.tag_surface == "unpaved":
        return False, "gravel"
    elif row.tag_surface in {"asphalt", "concrete"}:
        return True, row.tag_surface
    else:
        return True, row.tag_surface


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


def get_rehab_costs(
    row: pd.Series, rehab_costs: pd.DataFrame
) -> Tuple[float, float, str]:
    """
    Determine the cost of rehabilitation for a given road segment (row).

    Args:
        row (pd.Series): Road segment
        rehab_costs: (pd.DataFrame): Table of rehabilitation costs for various road types

    Returns:
        Tuple[float, float, str]: Minimum cost, maximum cost, units of cost
    """

    # bridge should be boolean after data cleaning step
    if row.bridge:
        highway_type = "bridge"
    else:
        highway_type = row.tag_highway

    if row.paved:
        condition = "paved"
    else:
        condition = "unpaved"

    minimum = rehab_costs.cost_min.loc[
        (rehab_costs.highway == highway_type) & (rehab_costs.road_cond == condition)
    ].squeeze()

    maximum = rehab_costs.cost_max.loc[
        (rehab_costs.highway == highway_type) & (rehab_costs.road_cond == condition)
    ].squeeze()

    unit = rehab_costs.cost_unit.loc[
        (rehab_costs.highway == highway_type) & (rehab_costs.road_cond == condition)
    ].squeeze()

    return float(minimum), float(maximum), str(unit)


def annotate_condition(
    network: snkit.network.Network, lane_width_m: float, shoulder_width_m: float
) -> snkit.network.Network:

    # infer paved status and material type from 'surface' column
    network.edges["paved_material"] = network.edges.apply(
        lambda x: get_road_condition(x), axis=1
    )
    # unpack 2 item iterable into two columns
    network.edges[["paved", "material"]] = network.edges["paved_material"].apply(
        pd.Series
    )

    # drop the now redundant columns
    network.edges.drop(["paved_material"], axis=1, inplace=True)

    # add number of lanes
    network.edges["lanes"] = network.edges.apply(lambda x: get_road_lanes(x), axis=1)

    # add road width
    network.edges["width_m"] = network.edges.apply(
        lambda x: x.lanes * lane_width_m + 2 * shoulder_width_m, axis=1
    )

    return network


def assign_road_speeds(row: pd.Series) -> Tuple[float, float]:
    """
    Infer road speed limits from road type and condition.

    Args:
        row (pd.Series): Road to infer speeds for, must have `highway`,
            `paved`, `Highway_min`, `Highway_max`, `Urban_min`, `Urban_max`,
            `Rural_min` and `Rural_max` labels.

    Returns:
        Tuple[float, float]: Likely minimum and maximum speeds
    """

    if row.tag_highway in {"motorway", "trunk", "primary"}:
        min_speed = float(row["Highway_min"])
        max_speed = float(row["Highway_max"])
    elif row.paved:
        min_speed = float(row["Urban_min"])
        max_speed = float(row["Urban_max"])
    else:
        min_speed = float(row["Rural_min"])
        max_speed = float(row["Rural_max"])

    # if we've got data from OSM, use that instead of the per-country data
    # use isnull because it works for numpy NaNs and Nones; a numpy NaN is truthy!
    if not pd.isnull(row.tag_maxspeed):
        max_speed = row.tag_maxspeed

    return min_speed, max_speed


def annotate_speeds(network: snkit.network.Network, speeds_by_country) -> snkit.network.Network:
    """
    Using OSM data (network.edges.maxspeed) and speeds_by_country, assemble a
    best guess of road speeds in km/h.

    Args:
        network (snkit.network.Network): Network to annotate. The network.edges
            must have a 'from_iso_a3' column
        speeds_by_country (pd.DataFrame): Speed limit information by country
            N.B. speeds_by_country is expected to have the following columns:
            'ISO_A3', 'CONTINENT', 'NAME', 'REGION', 'SUBREGION', 'Highway_min',
            'Highway_max', 'Rural_min', 'Rural_max', 'Urban_min', 'Urban_max'

    Returns:
        snkit.network.Network: Modified network
    """

    network.edges = pd.merge(
        network.edges,
        speeds_by_country,
        how="left",
        left_on=["from_iso_a3"],
        right_on=["ISO_A3"],
    )

    # infer a likely min and max road speed
    network.edges["min_max_speed"] = network.edges.apply(
        lambda x: assign_road_speeds(x), axis=1
    )

    # assign_road_speeds returned two values, unpack these into two columns
    network.edges[["min_speed_kmh", "max_speed_kmh"]] = network.edges[
        "min_max_speed"
    ].apply(pd.Series)

    # drop the intermediate columns
    network.edges.drop(
        ["min_max_speed"] + speeds_by_country.columns.values.tolist(),
        axis=1,
        inplace=True,
    )

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
            must have a 'from_iso_a3' column to merge on.
        transport_tariffs (pd.DataFrame): Table of transport tariffs by country and
            by mode of transport, with 'from_iso3', 'cost_km', 'cost_unit' and
            'cost_scaling' columns
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
        left_on=["from_iso_a3"],
        right_on=["from_iso3"],
    )
    network.edges["min_tariff"] = network.edges.apply(
        lambda x: float(x.tariff_cost) - (float(x.tariff_cost) * 0.2), axis=1
    )
    network.edges["max_tariff"] = network.edges.apply(
        lambda x: float(x.tariff_cost) + (float(x.tariff_cost) * 0.2), axis=1
    )
    network.edges.drop(["tariff_cost", "from_iso3"], axis=1, inplace=True)

    # calculate road segment lengths
    geod = Geod(ellps="WGS84")
    network.edges["length_m"] = network.edges.apply(
        lambda x: float(geod.geometry_length(x.geometry)), axis=1
    )

    # assign flow costs
    metres_per_km = 1_000

    network.edges["min_flow_cost"] = (
        flow_cost_time_factor * network.edges["length_m"] / metres_per_km
    ) / network.edges["max_speed_kmh"] + (
        network.edges["min_tariff"] * network.edges["length_m"] / metres_per_km
    )

    network.edges["max_flow_cost"] = (
        flow_cost_time_factor * network.edges["length_m"] / metres_per_km
    ) / network.edges["min_speed_kmh"] + (
        network.edges["max_tariff"] * network.edges["length_m"] / metres_per_km
    )

    network.edges["flow_cost_unit"] = "USD/ton"

    return network


if __name__ == "__main__":
    try:
        osm_edges_path = snakemake.input["edges"]  # type: ignore
        osm_nodes_path = snakemake.input["nodes"]  # type: ignore
        administrative_data_path = snakemake.input["admin"]  # type: ignore
        nodes_output_path = snakemake.output["nodes"]  # type: ignore
        edges_output_path = snakemake.output["edges"]  # type: ignore
        slice_number = snakemake.params["slice_number"]  # type: ignore
        road_speeds_path = snakemake.config["transport"]["speeds_path"]  # type: ignore
        rehabilitation_costs_path = snakemake.config["transport"]["rehabilitation_costs_path"]  # type: ignore
        transport_costs_path = snakemake.config["transport"]["tariff_costs_path"]  # type: ignore
        default_shoulder_width_metres = snakemake.config["transport"]["road"]["default_shoulder_width_metres"]  # type: ignore
        default_lane_width_metres = snakemake.config["transport"]["road"]["default_lane_width_metres"]  # type: ignore
        flow_cost_time_factor = snakemake.config["transport"]["road"]["flow_cost_time_factor"]  # type: ignore
        osm_epsg = snakemake.config["osm_epsg"]  # type: ignore
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
            road_speeds_path,
            rehabilitation_costs_path,
            transport_costs_path,
            default_shoulder_width_metres,
            default_lane_width_metres,
            flow_cost_time_factor,
            osm_epsg,
        ) = sys.argv[1:]
        # osm_edges_path = ../../results/geoparquet/tanzania-latest_filter-road/slice-0.geoparquet
        # osm_nodes_path = ../../results/geoparquet/tanzania-latest_filter-road/slice-0.geoparquet
        # administrative_data_path = ../../results/input/admin-boundaries/gadm36_levels.gpkg
        # nodes_output_path = ../../results/geoparquet/tanzania-latest_filter-road/slice-0_road_nodes.geoparquet
        # edges_output_path = ../../results/geoparquet/tanzania-latest_filter-road/slice-0_road_edges.geoparquet
        # slice_number = 0
        # road_speeds_path = ../../bundled_data/global_road_speeds.xlsx
        # rehabilitation_costs_path = ../../bundled_data/rehabilitation_costs.xlsx
        # transport_costs_path = ../../bundled_data/transport_costs.csv
        # default_shoulder_width_metres = 1.5
        # default_lane_width_metres = 3.25
        # flow_cost_time_factor = 0.49
        # osm_epsg = 4326

    # cast script arguments to relevant types where necessary
    default_shoulder_width_metres = float(default_shoulder_width_metres)
    default_lane_width_metres = float(default_lane_width_metres)
    flow_cost_time_factor = float(flow_cost_time_factor)
    slice_number = int(slice_number)
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
        write_empty_frames(edges_output_path, nodes_output_path)
        sys.exit(0)  # exit gracefully so snakemake will continue

    # osm_to_pq.py creates these columns but we're not using them, so discard
    edges = edges.drop(
        [col for col in edges.columns if col.startswith("start_node_") or col.startswith("end_node_")],
        axis="columns"
    )

    # for roads we do not currently use any nodes extracted from OSM (osm_nodes_path)
    logging.info("Creating road network")
    network = create_network(edges=clean_edges(edges), nodes=None, id_prefix=f"{slice_number}")
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
    network = annotate_condition(
        network, default_lane_width_metres, default_shoulder_width_metres
    )

    # select and label assets with their type
    # the asset_type is used to later select a damage curve
    network.edges.loc[network.edges.paved == False, 'asset_type'] = RoadAssets.UNPAVED
    network.edges.loc[network.edges.paved == True, 'asset_type'] = RoadAssets.PAVED
    network.edges.loc[network.edges.bridge == True, 'asset_type'] = RoadAssets.BRIDGE

    logging.info("Annotating network with road speed data")
    network = annotate_speeds(
        network, pd.read_excel(road_speeds_path, sheet_name="road")
    )

    logging.info("Annotating network with rehabilitation costs")
    network = annotate_rehabilitation_costs(
        network,
        pd.read_excel(rehabilitation_costs_path, sheet_name="road"),
        get_rehab_costs
    )

    logging.info("Annotating network with tariff and flow costs")
    network = annotate_tariff_flow_costs(
        network,
        pd.read_excel(transport_costs_path, sheet_name="road"),
        flow_cost_time_factor,
    )

    # TODO: drop superfluous columns? (e.g. OSM tags)

    # TODO: annotate with asset categories for direct damage estimation

    logging.info("Writing network to disk")
    network.edges.to_parquet(edges_output_path)
    network.nodes.to_parquet(nodes_output_path)

    logging.info("Done creating network")
