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


def get_administrative_data(file_path: str, to_epsg: int = None) -> gpd.GeoDataFrame:
    """
    Read administrative data (country ISO, country geometry) from disk

    Arguments:
        file_path (str): Location of file with country data
        to_epsg (int): EPSG code to project data to

    Returns:
        gpd.GeoDataFrame: Table of country and geometry data
    """

    # read file
    gdf = gpd.read_file(file_path)

    # check schema is as expected
    expected_columns = {'GID_0', 'NAME_0', 'geometry'}
    assert expected_columns == set(gdf.columns.values)

    # reproject if desired
    if to_epsg is not None:
        gdf = gdf.to_crs(epsg=to_epsg)

    # explode the ~250 country rows with multipolygon geometries into many thousands of rows with polygon geometries
    # ignore_index will forget the starting index and reindex the table
    gdf = gdf.explode(ignore_index=True)

    # rename these columns first so we don't have to do this twice (to nodes and edges) later
    gdf.rename(columns={"GID_0": "iso_code"}, inplace=True)

    # sort by and return
    return gdf.sort_values(by=["iso_code"], ascending=True)


def annotate_country(network: snkit.network.Network, countries: gpd.GeoDataFrame, crs_epsg: int) -> snkit.network.Network:
    """
    Label network edges and nodes with their country ISO code

    Arguments:
        network (snkit.network.Network): Network to label with geographic information
            network.edges should have 'edge_id', 'from_node_id', 'to_node_id', 'geometry'
            network.nodes should have 'node_id', 'geometry'
        countries (gpd.GeoDataFrame): Table expected to contain the following columns:
            'iso_code', 'NAME_0', 'geometry'
        crs_epsg (int): EPSG code for a standard CRS to use to compare geometries

    Returns
        snkit.network.Network: Labelled network
    """
    # unpack network for brevity in the following lines
    edges = network.edges
    nodes = network.nodes

    # CRS strings
    input_crs = f"EPSG:{crs_epsg}"
    projection_epsg = edges.estimate_utm_crs().to_epsg()
    projected_crs = f"EPSG:{projection_epsg}"
    logging.info(f"Inferred a suitable projection CRS of: {projected_crs}")

    # often we only want these columns
    core_node_columns = ["node_id", "iso_code", "geometry"]

    # spatial join nodes geometries to their containing country, retain only node geometries
    interior_nodes = gpd.sjoin(
        # subset nodes to only id
        nodes[["node_id", "geometry"]].to_crs(input_crs),
        countries.to_crs(input_crs),
        how="left",
        predicate='within'
    ).reset_index()
    interior_nodes = interior_nodes[~interior_nodes["iso_code"].isna()]
    interior_nodes = interior_nodes[core_node_columns]
    interior_nodes = interior_nodes.drop_duplicates(subset=["node_id"], keep="first")
    logging.info(f"Found {len(interior_nodes)} nodes that are within a country geometry")

    # for those nodes which weren't within a country geometry, label them with the nearest country
    # use a projected CRS so sjoin_nearest can more accurately compute distances
    # N.B. these projected UTM CRS are given for 6 degrees longitude sectors
    exterior_nodes = nodes[~nodes["node_id"].isin(interior_nodes["node_id"].values.tolist())]
    exterior_nodes = gpd.sjoin_nearest(
        exterior_nodes[["node_id", "geometry"]].to_crs(projected_crs),
        countries.to_crs(projected_crs),
        how="left"
    ).reset_index()
    exterior_nodes = exterior_nodes[core_node_columns]
    exterior_nodes = exterior_nodes.drop_duplicates(subset=["node_id"], keep="first")
    logging.info(
        f"Found {len(exterior_nodes)} nodes external to all country geometries, labelled "
        "these with their nearest country"
    )

    # join the outputs of the two joining processes together
    nodes = pd.concat([interior_nodes, exterior_nodes], axis=0, ignore_index=True)
    nodes = gpd.GeoDataFrame(
        nodes[core_node_columns],
        geometry="geometry",
        crs=input_crs
    )

    # set edge.from_node_id from node.node_id and use iso_code of from node as edge start
    edges = pd.merge(
        edges,
        nodes[["node_id", "iso_code"]],
        how="left",
        left_on=["from_node_id"],
        right_on=["node_id"]
    )
    edges.rename(columns={"iso_code": "from_iso"}, inplace=True)
    edges.drop("node_id", axis=1, inplace=True)

    # set edge.to_node_id from node.node_id and use iso_code of from node as edge end
    edges = pd.merge(
        edges,
        nodes[["node_id", "iso_code"]],
        how="left",
        left_on=["to_node_id"],
        right_on=["node_id"]
    )
    edges.rename(columns={"iso_code": "to_iso"}, inplace=True)
    edges.drop("node_id", axis=1, inplace=True)

    # encode country in id strings
    nodes["node_id"] = nodes.apply(lambda x: f"{x.iso_code}_{x.node_id}", axis=1)
    edges["from_node_id"] = edges.apply(lambda x: f"{x.from_iso}_{x.from_node_id}", axis=1)
    edges["to_node_id"] = edges.apply(lambda x: f"{x.to_iso}_{x.to_node_id}", axis=1)
    edges["edge_id"] = edges.apply(lambda x: f"{x.from_iso}_{x.to_iso}_{x.edge_id}", axis=1)

    network.nodes = nodes
    network.edges = edges

    return network


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

    if not row.surface:
        if row.highway in {'motorway', 'trunk', 'primary'}:
            return True, 'asphalt'
        else:
            return False, 'gravel'
    elif row.surface == 'paved':
        return True, 'asphalt'
    elif row.surface == 'unpaved':
        return False, 'gravel'
    elif row.surface in {'asphalt', 'concrete'}:
        return True, row.surface
    else:
        return True, row.surface


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
        lanes = int(row.lanes)
        if lanes < 1:
            return 1
        else:
            return lanes
    except (ValueError, TypeError):
        # couldn't cast row.lanes into an integer
        # instead guess at lanes from highway classification
        if row.highway in ('motorway', 'trunk', 'primary'):
            return 2
        else:
            return 1


def annotate_condition(
    network: snkit.network.Network,
    lane_width_m: float,
    shoulder_width_m: float
) -> snkit.network.Network:

    # calculate road segment lengths
    geod = Geod(ellps="WGS84")
    network.edges['length_m'] = network.edges.apply(
        lambda x: float(geod.geometry_length(x.geometry)),
        axis=1
    )

    # infer paved status and material type from 'surface' column
    network.edges['paved_material'] = network.edges.apply(
        lambda x: get_road_condition(x),
        axis=1
    )
    # unpack 2 item iterable into two columns
    network.edges[['paved', 'material']] = \
        network.edges['paved_material'].apply(pd.Series)

    # drop the now redundant columns
    network.edges.drop(['paved_material', 'surface'], axis=1, inplace=True)

    # add number of lanes
    network.edges['lanes'] = network.edges.apply(
        lambda x: get_road_lanes(x),
        axis=1
    )

    # add road width
    network.edges['width_m'] = network.edges.apply(
        lambda x: x.lanes * lane_width_m + 2 * shoulder_width_m,
        axis=1
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

    if row.highway in {'motorway', 'trunk', 'primary'}:
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
    if not pd.isnull(row.maxspeed):
        max_speed = row.maxspeed

    return min_speed, max_speed


def annotate_speeds(
    network: snkit.network.Network,
    speeds_by_country
) -> snkit.network.Network:
    """
    Using OSM data (network.edges.maxspeed) and speeds_by_country, assemble a
    best guess of road speeds in km/h.

    Args:
        network (snkit.network.Network): Network to annotate. The network.edges
            must have a 'from_iso' column
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
        left_on=["from_iso"],
        right_on=["ISO_A3"]
    )

    # infer a likely min and max road speed
    network.edges["min_max_speed"] = network.edges.apply(
        lambda x: assign_road_speeds(x),
        axis=1
    )

    # assign_road_speeds returned two values, unpack these into two columns
    network.edges[["min_speed_kmh", "max_speed_kmh"]] = \
        network.edges["min_max_speed"].apply(pd.Series)

    # drop the intermediate columns
    network.edges.drop(
        ["min_max_speed", "maxspeed"] + speeds_by_country.columns.values.tolist(),
        axis=1,
        inplace=True
    )

    return network


def get_rehab_costs(row: pd.Series, rehab_costs: pd.DataFrame) -> Tuple[float, float, str]:
    """
    Determine the cost of rehabilitation for a given road segment (row).

    Args:
        row (pd.Series): Road segment
        rehab_costs: (pd.DataFrame): Table of rehabilitation costs for various road types

    Returns:
        Tuple[float, float, str]: Minimum cost, maximum cost, units of cost
    """

    # bridge should be a boolean type after data cleaning step
    if row.bridge:
        highway_type = "bridge"
    else:
        highway_type = row.highway

    if row.paved:
        condition = 'paved'
    else:
        condition = 'unpaved'

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


def annotate_rehabilitation_costs(
    network: snkit.network.Network,
    rehab_costs: pd.DataFrame
) -> snkit.network.Network:

    # lookup costs
    network.edges["rehab_costs"] = network.edges.apply(
        get_rehab_costs,
        axis=1,
        args=(rehab_costs,)
    )

    # unpack results into 3 columns
    network.edges[["rehab_cost_min", "rehab_cost_max", "rehab_cost_unit"]] = \
        network.edges["rehab_costs"].apply(pd.Series)
    network.edges.drop(["rehab_costs"], axis=1, inplace=True)

    return network


def annotate_tariff_flow_costs(
    network: snkit.network.Network,
    transport_tariffs: pd.DataFrame,
    transport_type: str,
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
        transport_type (str): Transport category, one of 'road', 'rail' or 'IWW'
            (internal waterways).
        flow_cost_time_factor (float): A fudge factor that varies (by country?)
            This may well need consuming as location specific data in future.
    Returns:
        snkit.network.Network: Modified network
    """

    # input checking
    expected_columns = {
        'from_iso3', 'transport', 'cost_km', 'cost_unit', 'cost_scaling'
    }
    if not set(transport_tariffs.columns).issuperset(expected_columns):
        raise ValueError(f"{expected_columns=} for transport_tariffs")
    expected_transport_types = {'road', 'rail', 'IWW'}
    if transport_type not in expected_transport_types:
        raise ValueError(f"{transport_type=} not in {expected_transport_types=}")

    # rename and subset table
    transport_tariffs.rename(columns={'cost_km': 'tariff_cost', 'cost_unit': 'tariff_unit'}, inplace=True)
    transport_tariffs = transport_tariffs.loc[transport_tariffs['transport'] == "road"]

    # merge datasets
    network.edges = pd.merge(
        network.edges,
        transport_tariffs[["from_iso3", "tariff_cost", "tariff_unit"]],
        how="left",
        left_on=["from_iso"],
        right_on=["from_iso3"]
    )
    network.edges["min_tariff"] = network.edges.apply(
        lambda x: float(x.tariff_cost) - (float(x.tariff_cost) * 0.2),
        axis=1
    )
    network.edges["max_tariff"] = network.edges.apply(
        lambda x: float(x.tariff_cost) + (float(x.tariff_cost) * 0.2),
        axis=1
    )
    network.edges.drop(["tariff_cost", "from_iso3"], axis=1, inplace=True)

    # assign flow costs
    metres_per_km = 1_000

    network.edges["min_flow_cost"] = \
        (flow_cost_time_factor * network.edges["length_m"] / metres_per_km) \
        / network.edges["max_speed_kmh"] \
        + (network.edges["min_tariff"] * network.edges["length_m"] / metres_per_km)

    network.edges["max_flow_cost"] = \
        (flow_cost_time_factor * network.edges["length_m"] / metres_per_km) \
        / network.edges["min_speed_kmh"] \
        + (network.edges["max_tariff"] * network.edges["length_m"] / metres_per_km)

    network.edges["flow_cost_unit"] = 'USD/ton'

    return network


if __name__ == '__main__':
    try:
        nodes_path = snakemake.input["nodes"]
        edges_path = snakemake.input["edges"]
        administrative_data_path = snakemake.input["admin"]
        output_nodes_path = snakemake.output["nodes"]
        output_edges_path = snakemake.output["edges"]
        road_speeds_path = snakemake.config["road_speeds_path"]
        rehabilitation_costs_path = snakemake.config["road_rehabilitation_costs_path"]
        transport_costs_path = snakemake.config["transport_costs_path"]
        default_shoulder_width_metres = snakemake.config["road_default_shoulder_width_metres"]
        default_lane_width_metres = snakemake.config["road_default_lane_width_metres"]
        flow_cost_time_factor = snakemake.config["road_flow_cost_time_factor"]
        osm_epsg = snakemake.config["osm_epsg"]
    except NameError:
        # If "snakemake" doesn't exist then must be running from the
        # command line.
        nodes_path, edges_path, administrative_data_path, output_nodes_path, output_edges_path, \
            road_speeds_path, rehabilitation_costs_path, transport_costs_path, \
            default_shoulder_width_metres, default_lane_width_metres, flow_cost_time_factor, \
            osm_epsg = sys.argv[1:]
        # nodes_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_nodes.geoparquet
        # edges_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_edges.geoparquet
        # administrative_data_path = ../../results/input/admin-boundaries/gadm36_levels.gpkg
        # output_nodes_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_nodes_annotated.geoparquet
        # output_edges_path = ../../results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_edges_annotated.geoparquet
        # road_speeds_path = ../../bundled_data/global_road_speeds.xlsx
        # rehabilitation_costs_path = ../../bundled_data/rehabilitation_costs.xlsx
        # transport_costs_path = ../../bundled_data/transport_costs.csv
        # default_shoulder_width_metres = 1.5
        # default_lane_width_metres = 6.5
        # flow_cost_time_factor = 0.49
        # osm_epsg = 4326

    # cast script arguments to numeric types where necessary
    default_shoulder_width_metres = float(default_shoulder_width_metres)
    default_lane_width_metres = float(default_lane_width_metres)
    flow_cost_time_factor = float(flow_cost_time_factor)
    osm_epsg = int(osm_epsg)

    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

    logging.info('Reading network from disk')
    warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')
    try:
        edges = gpd.read_parquet(edges_path)
        nodes = gpd.read_parquet(nodes_path)
    except ValueError as error:
        # if the input parquet file does not contain a geometry column, geopandas
        # will raise a ValueError rather than try to procede
        logging.info(
            f"{error}\n"
            "writing empty files and skipping processing..."
        )

        # snakemake requires that output files exist though, so write empty ones
        empty_gdf = gpd.GeoDataFrame([])
        empty_gdf.to_parquet(output_edges_path)
        empty_gdf.to_parquet(output_nodes_path)
        sys.exit(0)  # exit gracefully so snakemake will continue

    annotated_network = snkit.network.Network(
        edges=gpd.read_parquet(edges_path),
        nodes=gpd.read_parquet(nodes_path),
    )
    logging.info(
        f'Network contains {len(annotated_network.edges)} edges and {len(annotated_network.nodes)} nodes'
    )

    logging.info('Annotating network with administrative data')
    annotated_network = annotate_country(
        annotated_network,
        get_administrative_data(administrative_data_path, to_epsg=osm_epsg),
        osm_epsg
    )

    logging.info('Annotating network with road type and condition data, calculated segment lengths')
    annotated_network = annotate_condition(
        annotated_network,
        default_lane_width_metres,
        default_shoulder_width_metres
    )

    logging.info('Annotating network with road speed data')
    annotated_network = annotate_speeds(
        annotated_network,
        pd.read_excel(road_speeds_path, sheet_name="global speeds")
    )

    logging.info('Annotating network with rehabilitation costs')
    annotated_network = annotate_rehabilitation_costs(
        annotated_network,
        pd.read_excel(rehabilitation_costs_path, sheet_name="road_costs")
    )

    logging.info('Annotating network with tariff and flow costs')
    annotated_network = annotate_tariff_flow_costs(
        annotated_network,
        pd.read_csv(transport_costs_path),
        "road",
        flow_cost_time_factor,
    )

    logging.info('Writing network to disk')
    annotated_network.edges.to_parquet(output_edges_path)
    annotated_network.nodes.to_parquet(output_nodes_path)

    logging.info('Done annotating network with flow and rehabilitation data')
