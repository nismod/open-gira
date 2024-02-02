"""
For each storm (maximum wind speed field) in storm list:
- Fail edge segments who experience a wind speed greater than given threshold
- Record length of edges in exceedence of threshold
- Attempt to allocate power from sources to sinks over degraded network
- Calculate ratio of nominal power to degraded power, the 'supply factor'
- Estimate number of customers affected
"""

import json
import logging
import os
import multiprocessing
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

from open_gira import fields
from open_gira.grid import weighted_allocation
from open_gira.io import bit_pack_dataset_encoding, bit_pack_dataarray_encoding
import snkit


# do not store supply_factor values greater than this
# a value less than 1 limits the targets stored in disruption files to areas
# negatively affected by storms
MAX_SUPPLY_FACTOR: float = 0.90


def process_event(
    storm_id: str,
    wind_fields: xr.DataArray,
    splits: gpd.GeoDataFrame,
    speed_thresholds: list[float],
    network: snkit.network.Network,
    exposure_dir: str,
    disruption_dir: str
) -> None:
    """
    File handling, filtering, logging and serialisation to disk for the exposure
    and disruption associated with a given storm.

    Args:
        storm_id: Unique ID of event to process
        wind_fields: Max wind speeds on event_id, longitude, latitude coordinates
        splits: Infrastructure edges split by wind grid
        speed_thresholds: Wind speeds at which to consider assets failed
        network: Network of nodes and edges to allocate power over
        exposure_dir: Location to write out exposure (length of assets exposed
            to wind speeds in excess of threshold(s)).
        disruption_dir: Location to write out disruption (supply factors and
            number of people disconnected)
    """

    logging.info(storm_id)
    exposure, disruption = degrade_grid_with_storm(storm_id, wind_fields, splits, speed_thresholds, network)

    # filter out values (and coordinate values) that are not of interest
    # this helps keep output file sizes manageable (especially for large networks)
    exposure = exposure.sel(event_id=storm_id)
    exposure = exposure.where(exposure.length_m > 0, drop=True)
    disruption = disruption.sel(event_id=storm_id)
    disruption = disruption.where(disruption.supply_factor < MAX_SUPPLY_FACTOR, drop=True)

    exposure_summary = exposure.length_m.sum(dim='edge')
    exposure_summary_str = (
        "Exposure summary" +
        "\nThreshold [m s-1], Grid exposed [m]\n" +
        "\n".join([f"{exposure.threshold:.1f}, {exposure:.2E}" for exposure in exposure_summary])
    )
    logging.info(exposure_summary_str)

    logging.info("Writing results to disk")

    # pack floating point data into 16 bit integer types to save space on disk
    exposure.to_netcdf(
        os.path.join(exposure_dir, f"{storm_id}.nc"),
        encoding=bit_pack_dataset_encoding(exposure)
    )

    # supply factor is a float between 0 and 1, can reasonably pack into 16 bit integer space
    disruption_encoding = bit_pack_dataarray_encoding(disruption.supply_factor)

    # customers affected spans a large (positive) range, use a single precision
    # float to save space over the default double
    disruption_encoding["customers_affected"] = {"dtype": "float32"}

    disruption.to_netcdf(
        os.path.join(disruption_dir, f"{storm_id}.nc"),
        encoding=disruption_encoding
    )

    return


def subset_network(
    nominal_network: snkit.network.Network,
    failed_edge_ids: np.ndarray,
    surviving_edge_ids: np.ndarray,
    id_col: str = "component_id"
) -> snkit.network.Network:
    """
    Take `nominal_network`, remove failed edges. Relabel nodes and edges of the
    resulting network with component ids, in `id_col` column.

    Assumes:
        - edge IDs are edges' index in the network.edge table
        - node IDs are nodes' index in the network.node table

    Arguments:
        nominal_network: Intact network
        failed_edge_ids: Edges to exclude from a new network
        surviving_edge_ids: Edges with which to create a new network
        id_col: Name of column to store component ids

    Returns:
        Degraded network, labelled with component ids
    """

    n_edges_nominal = len(failed_edge_ids) + len(surviving_edge_ids)
    assert n_edges_nominal == len(nominal_network.edges)

    # construct network from what remains
    degraded_network = snkit.network.Network(
        edges=nominal_network.edges.loc[surviving_edge_ids].copy(),
        nodes=nominal_network.nodes.copy()
    )

    connected_components: list[set] = list(snkit.network.get_connected_components(degraded_network))

    edge_component_ids: np.ndarray = np.zeros(n_edges_nominal, dtype=int)
    node_component_ids: np.ndarray = np.zeros(len(nominal_network.nodes), dtype=int)

    # record which edges we've already checked, so we can reduce the size of our `isin` query
    checked_edges = np.full(n_edges_nominal, False)

    # don't want to check any edge that's not in the surviving network
    checked_edges[failed_edge_ids] = True

    for count, component in enumerate(connected_components):

        component_id: int = count + 1

        # boolean series of component membership
        component_edge_mask: pd.Series = nominal_network.edges.loc[~checked_edges, "from_id"].isin(component)

        # update our record of which edges we've checked
        checked_edges[component_edge_mask.index.values] = component_edge_mask.values

        # indicies of all edges in this component
        component_member_indicies: np.ndarray = component_edge_mask[component_edge_mask].index.values

        # assign the component id
        edge_component_ids[component_member_indicies] = component_id

        # index into node component ids using the ids from the component, set to component id
        node_component_ids[np.fromiter(component, int, len(component))] = component_id

    # assign results
    degraded_network.edges[id_col] = edge_component_ids[surviving_edge_ids]
    degraded_network.nodes[id_col] = node_component_ids

    assert all(degraded_network.edges[id_col] != 0)
    assert all(degraded_network.nodes[id_col] != 0)

    return degraded_network


def build_dataset(var_names: tuple[str], dim_type: dict[str, type], **kwargs) -> xr.Dataset:
    """
    Build an empty (NaN filled) xarray Dataset given names, types and coordinate values.

    Args:
        dim_type: Coordinate dimension name to type mapping, e.g. {"event_id": str, "target": int}
        var_names: Names of variables to create across all dimensions. Filled with np.nan.

        Additional kwargs:
            Each key in `dim_type` must be matched by an identically named
            kwarg pointing to the iterable to use for the coordinate dimension
            values. For the example mapping, expect to be called as follows:

            build_dataset(
                {"event_id": str, "target": int},
                event_id=["2018N123"],
                target=[1, 2, 3]
            )

            If no matching kwarg is found, an empty coordinate array will be used.

    Returns:
        Dataset with coordinates dimensions as specified, with named variables filled with np.nan.
    """

    dim_coords = {}
    for dim_name in dim_type.keys():
        try:
            dim_coords[dim_name] = kwargs[dim_name]
        except KeyError:
            # if we haven't been passed coords, use an empty list
            dim_coords[dim_name] = []

    return xr.Dataset(
        data_vars={
            var_name:
            (
                dim_type.keys(),
                np.full(tuple(len(coord) for coord in dim_coords.values()), np.nan)
            )
            for var_name in var_names
        },
        coords={
            dim_name: np.array(np.atleast_1d(dim_coords[dim_name]), dtype=dim_type[dim_name])
            for dim_name in dim_type
        }
    )


def degrade_grid_with_storm(
    storm_id: str,
    wind_fields: xr.DataArray,
    splits: gpd.GeoDataFrame,
    speed_thresholds: list,
    network: snkit.network.Network
) -> tuple[xr.Dataset, xr.Dataset]:
    """
    Use a maximum wind speed field and a electricity grid representation,
    degrade the network for a set of damage speed thresholds. Estimate the
    reduction in available supply and number of customers affected.

    Args:
        storm_id: String ID of storm to simulate
        wind_fields: Maximum wind speeds experienced in gridded domain.
        splits: Electricity grid split over raster grid. Note that `splits`
            should contain two columns to positionally index `wind_fields`
        speed_thresholds: List of wind speeds to fail network edges at. Should
            be in the same units as `wind_fields`.
        network: Network representation of electricity grid. Edges should have
            topology. Nodes should have an `asset_type` and where `asset_type`
            is 'target', there should be a nominal power consumption allocated.

    Returns:
        Dataset containing length_m exposure variable on event_id, threshold
            and edge dimensions.
        Dataset containing supply_factor and customers_affected disruption
            variables on event_id, threshold and target dimensions.
    """

    # N.B. we have a generic node 'id' but also a 'target_id' which should only
    # be set for target nodes -- it is globally unique and corresponds to the
    # global targets file (which contains their vector geometry)
    try:
        target_ids = network.nodes[network.nodes.asset_type == "target"].target_id.astype(int).values
    except AttributeError:
        logging.info("No viable network available, returning null result.")
        return (
            build_dataset(("length_m",), {"event_id": str, "threshold": float, "edge": int}, event_id=[storm_id]),
            build_dataset(("supply_factor", "customers_affected"), {"event_id": str, "threshold": float, "target": int}, event_id=[storm_id])
        )

    exposure = build_dataset(
        ("length_m",),
        {"event_id": str, "threshold": float, "edge": int},
        event_id=[storm_id],
        threshold=speed_thresholds,
        edge=network.edges.id.astype(int).values
    )

    disruption = build_dataset(
        ("supply_factor", "customers_affected"),
        {"event_id": str, "threshold": float, "target": int},
        event_id=[storm_id],
        threshold=speed_thresholds,
        target=target_ids
    )

    try:
        # rank 1, length of splits DataFrame
        # N.B. to index at points rather than the cross-product of indicies, index with DataArrays
        # https://docs.xarray.dev/en/stable/user-guide/indexing.html#vectorized-indexing
        max_wind_speeds: xr.DataArray = wind_fields.sel(event_id=storm_id).isel(
            longitude=splits[fields.RASTER_I].to_xarray(),
            latitude=splits[fields.RASTER_J].to_xarray()
        )
    except KeyError:
        logging.info("No wind field available, returning null result.")
        return exposure, disruption

    # index and `id` column need to match as we will select rows by indexing with ids
    assert all(network.edges.index == network.edges.id)

    # sort into ascending order; if no damage at a given threshold,
    # more resilient thresholds are guaranteed to be safe
    for threshold in sorted(speed_thresholds):
        survival_mask: pd.Series = (max_wind_speeds < threshold).to_pandas().loc[:, "max_wind_speed"]

        try:
            n_failed: int = survival_mask.value_counts()[False]
        except KeyError:
            # there is no damage, return early
            logging.info(f"No damage detected at {threshold} ms-1")
            return exposure, disruption

        ############
        # EXPOSURE #
        ############

        # all splits above threshold
        failed_splits: gpd.GeoDataFrame = splits.set_index("id", drop=True).loc[~survival_mask].copy()
        if failed_splits.empty is True:
            return exposure, disruption
        # label failed splits with length
        failed_splits["length_m"] = failed_splits.to_crs(failed_splits.estimate_utm_crs()).geometry.length
        # sum across edge id to find exposed length in case where line split
        # reset_index gives us our edge id column back
        exposed_edge_lengths = failed_splits[["length_m"]].groupby("id").sum().reset_index()
        # store result in dataset
        indicies = dict(event_id=storm_id, threshold=threshold, edge=exposed_edge_lengths.id.astype(int).values)
        exposure.length_m.loc[indicies] = exposed_edge_lengths.length_m

        ##############
        # DISRUPTION #
        ##############

        # edge ids containing a split below wind speed threshold
        failed_edge_ids: np.ndarray = splits.loc[~survival_mask, "id"].unique()
        surviving_edge_ids: np.ndarray = network.edges.loc[~network.edges.id.isin(failed_edge_ids), "id"].values

        # construct network from what remains, relabel connectedness
        surviving_network: snkit.network.Network = subset_network(network, failed_edge_ids, surviving_edge_ids)
        c_nominal: int = len(set(network.nodes.component_id))
        c_surviving: int = len(set(surviving_network.nodes.component_id))

        fraction_failed: float = n_failed / len(survival_mask)
        failure_str = "{:s} -> {:.2f}% edges failed @ {:.1f} [m/s] threshold, {:d} -> {:d} components"
        logging.info(failure_str.format(str(storm_id), 100 * fraction_failed, threshold, c_nominal, c_surviving))

        # about to mutate power_mw column, make a copy first
        surviving_network.nodes["power_nominal_mw"] = surviving_network.nodes["power_mw"]

        # if there's no gdp data available at all, use the population as a weight
        # this should have be used when creating the network in create_electricity_network.py
        if surviving_network.nodes[surviving_network.nodes.asset_type == "target"].gdp.sum() == 0:
            weight_col = "population"
        else:
            weight_col = "gdp"

        # allocate power within components, from sources to targets, weighted (typically) by gdp of targets
        targets: pd.DataFrame = weighted_allocation(
            surviving_network.nodes,
            variable_col="power_mw",
            weight_col=weight_col,
            component_col="component_id",
            asset_col="asset_type",
            source_name="source",
            sink_name="target",
        )

        # calculate ratio of degraded power supply to nominal power supply
        targets["supply_factor"] = targets.power_mw / targets.power_nominal_mw

        # calculate the number of customers affect in each target
        # N.B. supply_factor can be > 1
        # so clip to zero to prevent negative customers_affected in areas with 'oversupply'
        targets["customers_affected"] = np.clip(1 - targets.supply_factor, 0, None) * targets.population

        # assign data to dataset
        indicies = dict(event_id=storm_id, threshold=threshold, target=targets.target_id.astype(int).values)
        disruption.supply_factor.loc[indicies] = targets.supply_factor
        disruption.customers_affected.loc[indicies] = targets.customers_affected

    return exposure, disruption


if __name__ == "__main__":

    edges_path: str = snakemake.input.grid_edges
    nodes_path: str = snakemake.input.grid_nodes
    splits_path: str = snakemake.input.grid_splits
    wind_speeds_path: str = snakemake.input.wind_speeds
    speed_thresholds: list[float] = snakemake.config["transmission_windspeed_failure"]
    exposure_dir: str = snakemake.output.exposure
    disruption_dir: str = snakemake.output.disruption
    n_proc: int = int(snakemake.threads)

    os.makedirs(exposure_dir)
    os.makedirs(disruption_dir)

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    logging.info("Loading wind speed metadata")
    wind_fields: xr.Dataset = xr.open_dataset(wind_speeds_path)

    if len(wind_fields.variables) == 0:
        logging.info("Empty wind speed file, skipping...")
        sys.exit(0)

    logging.info(wind_fields.max_wind_speed)

    logging.info("Loading network data")
    network = snkit.network.Network(
        edges=gpd.read_parquet(edges_path),
        nodes=gpd.read_parquet(nodes_path)
    )
    splits: gpd.GeoDataFrame = gpd.read_parquet(splits_path).set_crs(epsg=4326)
    logging.info(f"{len(network.edges)} network edges")
    logging.info(f"{len(network.nodes)} network nodes")

    logging.info(f"Using damage thresholds: {speed_thresholds} [m s-1]")

    args = (
        (storm_id.item(), wind_fields, splits, speed_thresholds, network, exposure_dir, disruption_dir)
        for storm_id in wind_fields.event_id
    )
    if n_proc > 1:
        with multiprocessing.get_context("fork").Pool(processes=n_proc) as pool:
            pool.starmap(process_event, args)
    else:
        for arg in args:
            process_event(*arg)
