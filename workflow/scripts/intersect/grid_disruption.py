"""
For a given storm (maximum wind speed field):
    1) Fail edge segments who experience a wind speed greater than given threshold
    2) Attempt to allocate power from sources to sinks over degraded network
    3) Calculate ratio of nominal power to degraded power, the 'supply factor'
    4) Estimate number of customers affected
"""

import logging
from typing import Optional
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

from open_gira import fields
from open_gira.grid import weighted_allocation
import snkit


# coordinates that the exposure variables can be indexed by
EXPOSURE_COORDS: dict[str, type] = {
    "event_id": str,
    "threshold": float,
    "target": int,
}

NETCDF_ENCODING = {variable: {"zlib": True, "complevel": 9} for variable in EXPOSURE_COORDS.keys()}


def exposure_dataset(
    event_id: Optional[list[str]] = None,
    thresholds: Optional[list[float]] = None,
    targets: Optional[list[int]] = None,
) -> xr.Dataset:
    """
    Build results object from given coordinate arguments.

    Args:
        event_id: Storm identifier
        thresholds: Iterable of wind speed damage thresholds
        targets: Iterable of globally unique target IDs
    """

    # if we haven't been passed coords, use empty lists
    if event_id is None:
        event_id = []
    if thresholds is None:
        thresholds = []
    if targets is None:
        targets = []

    shape = (len(event_id), len(thresholds), len(targets))
    return xr.Dataset(
        data_vars=dict(
            supply_factor=(EXPOSURE_COORDS.keys(), np.full(shape, np.nan)),
            customers_affected=(EXPOSURE_COORDS.keys(), np.full(shape, np.nan))
        ),
        coords=dict(
            # scalar dimensions result in ValueError, use atleast_1d as workaround
            # https://stackoverflow.com/a/58858160
            event_id=np.atleast_1d(event_id),
            threshold=np.atleast_1d(thresholds),
            target=np.atleast_1d(targets)
        )
    )


def degrade_grid_with_storm(
    storm_id: xr.DataArray,
    wind_fields: xr.DataArray,
    splits: pd.DataFrame,
    speed_thresholds: list,
    network: snkit.network.Network
) -> xr.Dataset:
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
        Dataset containing supply_factor and customers_affected variables on
            event_id, threshold and target dimensions.
    """

    # N.B. we have a generic node 'id' but also a 'target_id' which should only
    # be set for target nodes -- it is globally unique and corresponds to the
    # global targets file (which contains their vector geometry)
    try:
        target_ids = network.nodes[network.nodes.asset_type == "target"].target_id.astype(int).values
    except AttributeError:
        logging.info("No viable network available, returning null result.")
        return exposure_dataset(event_id=[storm_id], thresholds=speed_thresholds)

    # build coordinates for results object
    exposure = exposure_dataset(event_id=[storm_id], thresholds=speed_thresholds, targets=target_ids)

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
        return exposure

    # sort into ascending order; if no damage at a given threshold,
    # more resilient thresholds are guaranteed to be safe
    for threshold in sorted(speed_thresholds):
        survival_mask: pd.Series = (max_wind_speeds < threshold).to_pandas()

        try:
            n_failed: int = survival_mask.value_counts()[False]
        except KeyError:
            # there is no damage, return early
            return exposure

        surviving_edge_ids = set(splits.loc[survival_mask, "id"])
        surviving_edges: pd.DataFrame = network.edges.loc[network.edges.id.isin(surviving_edge_ids), :]
        surviving_network = snkit.network.Network(
            edges=surviving_edges.copy(),
            nodes=network.nodes.copy()
        )

        # check topology of degraded network
        surviving_network = snkit.network.add_component_ids(surviving_network)
        c_nominal = len(set(network.nodes.component_id))
        c_surviving = len(set(surviving_network.nodes.component_id))

        fraction_failed: float = n_failed / len(survival_mask)
        failure_str = "{:s} -> {:.2f}% edges failed @ {:.1f} [m/s] threshold, {:d} -> {:d} components"
        logging.info(failure_str.format(str(storm_id.values), 100 * fraction_failed, threshold, c_nominal, c_surviving))

        # about to mutate power_mw column, make a copy first
        surviving_network.nodes["power_nominal_mw"] = surviving_network.nodes["power_mw"]

        # allocate power within components, from sources to targets, weighted by gdp of targets
        targets: pd.DataFrame = weighted_allocation(
            surviving_network.nodes,
            variable_col="power_mw",
            weight_col="gdp",
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
        exposure.supply_factor.loc[indicies] = targets.supply_factor
        exposure.customers_affected.loc[indicies] = targets.customers_affected

    return exposure


if __name__ == "__main__":

    edges_path: str = snakemake.input.grid_edges
    nodes_path: str = snakemake.input.grid_nodes
    splits_path: str = snakemake.input.grid_splits
    wind_speeds_path: str = snakemake.input.wind_speeds
    speed_thresholds: list[float] = snakemake.config["transmission_windspeed_failure"]
    storm_id: str = snakemake.wildcards.STORM_ID
    exposure_path: str = snakemake.output.exposure

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    logging.info("Loading wind speed data")
    wind_fields: xr.Dataset = xr.open_dataset(wind_speeds_path)
    if len(wind_fields.variables) == 0:
        logging.info("Empty wind speed file, writing null result to disk.")
        exposure = exposure_dataset(event_id=[storm_id], thresholds=speed_thresholds)
        exposure.to_netcdf(exposure_path, encoding=NETCDF_ENCODING)
        sys.exit(0)

    logging.info(wind_fields.max_wind_speed)  # use xarray repr

    logging.info("Loading network data")
    network = snkit.network.Network(
        edges=gpd.read_parquet(edges_path),
        nodes=gpd.read_parquet(nodes_path)
    )
    splits: pd.DataFrame = pd.read_parquet(splits_path)
    logging.info(f"{len(network.edges)} network edges")
    logging.info(f"{len(network.nodes)} network nodes")

    logging.info(f"Using damage thresholds: {speed_thresholds} [m/s]")

    logging.info("Simulating electricity network failure due to wind damage...")
    exposure = degrade_grid_with_storm(storm_id, wind_fields, splits, speed_thresholds, network)

    logging.info(f"Writing results to disk")
    exposure.to_netcdf(exposure_path, encoding=NETCDF_ENCODING)
