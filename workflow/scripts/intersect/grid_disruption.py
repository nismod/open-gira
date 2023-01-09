"""
For each maximum wind field associated with a storm:
    1) Fail edge segments who experience a wind speed greater than given threshold
    2) Attempt to allocate power from sources to sinks over degraded network
    3) Calculate ratio of nominal power to degraded power, the 'supply factor'
    4) Estimate number of customers affected
"""

import logging
import multiprocessing

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

from open_gira import fields
from open_gira.grid import allocate_power_to_targets
import snkit


logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

DIM_NAMES = ["event_id", "threshold", "target"]


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

    storm_id_str = str(storm_id.values)  # coerce singleton array into str
    logging.info(storm_id_str)

    # rank 1, length of splits DataFrame
    # N.B. to index at points rather than the cross-product of indicies, index with DataArrays
    # https://docs.xarray.dev/en/stable/user-guide/indexing.html#vectorized-indexing
    max_wind_speeds: xr.DataArray = wind_fields.sel(event_id=storm_id).isel(
        long=splits[fields.RASTER_I].to_xarray(),
        lat=splits[fields.RASTER_J].to_xarray()
    )

    # object for collating results from each damage threshold
    target_ids = network.nodes[network.nodes.asset_type == "target"].id.values
    return_shape = (1, len(speed_thresholds), len(target_ids))
    empty = np.full(return_shape, np.nan)
    exposure = xr.Dataset(
        data_vars=dict(
            supply_factor=(DIM_NAMES, empty),
            customers_affected=(DIM_NAMES, empty)
        ),
        coords=dict(
            # scalar dimensions result in ValueError, use atleast_1d as workaround
            # https://stackoverflow.com/a/58858160
            event_id=np.atleast_1d(storm_id),
            threshold=speed_thresholds,
            target=target_ids
        )
    )

    for threshold in speed_thresholds:
        survival_mask: pd.Series = (max_wind_speeds < threshold).to_pandas()

        try:
            n_failed: int = survival_mask.value_counts()[False]
        except KeyError:
            # network is intact, short circuit
            return exposure

        surviving_edge_ids = set(splits.loc[survival_mask, "id"])
        surviving_edges: pd.DataFrame = network.edges.loc[network.edges.id.isin(surviving_edge_ids), :]
        surviving_network = snkit.network.Network(
            edges=surviving_edges.copy(),
            nodes=network.nodes.copy()
        )

        # check topology of degraded network
        surviving_network = snkit.network.add_topology(surviving_network, id_col="id")
        surviving_network = snkit.network.add_component_ids(surviving_network)
        c_nominal = len(set(network.nodes.component_id))
        c_surviving = len(set(surviving_network.nodes.component_id))

        fraction_failed: float = n_failed / len(survival_mask)
        failure_str = "{:s} -> {:>6.2f}% edges failed @ {:.1f} [m/s] threshold, {:d} -> {:d} components"
        logging.info(failure_str.format(storm_id_str, 100 * fraction_failed, threshold, c_nominal, c_surviving))

        # reallocate generating capacity by GDP within components
        nodes = surviving_network.nodes
        # about to mutate power_mw column, make a copy first
        nodes["power_nominal_mw"] = nodes["power_mw"]
        # for new network configuration, allocate source generation capacity to targets, weighted by GDP
        nodes = allocate_power_to_targets(nodes, "gdp")

        # calculate ratio of degraded power supply to nominal power supply
        targets = nodes[nodes.asset_type == "target"].copy()
        targets["supply_factor"] = targets.power_mw / targets.power_nominal_mw

        # supply factor can be > 1, so clip to zero to prevent negative customers_affected in areas with 'oversupply'
        targets["customers_affected"] = np.clip(1 - targets.supply_factor, 0, None) * targets.population

        # assign data to dataset
        indicies = dict(event_id=storm_id, threshold=threshold, target=targets.id.values)

        exposure.supply_factor.loc[indicies] = targets.supply_factor
        exposure.customers_affected.loc[indicies] = targets.customers_affected

    return exposure


if __name__ == "__main__":

    edges_path: str = snakemake.input.grid_edges
    nodes_path: str = snakemake.input.grid_nodes
    splits_path: str = snakemake.input.grid_splits
    wind_speeds_path: str = snakemake.input.wind_speeds
    # sort into ascending order; if no damage at a given threshold,
    # more resilient thresholds are guaranteed to be safe
    speed_thresholds: list[float] = sorted(snakemake.config["transmission_windspeed_failure"])
    parallel: bool = snakemake.config["parallelise_by_storm"]
    damages_path: str = snakemake.output.damages

    logging.info("Loading network data")

    network = snkit.network.Network(
        edges=gpd.read_parquet(edges_path),
        nodes=gpd.read_parquet(nodes_path)
    )
    logging.info(f"{len(network.edges)} network edges")
    logging.info(f"{len(network.nodes)} network nodes")

    splits: pd.DataFrame = pd.read_parquet(splits_path)

    logging.info("Loading wind speed data")

    wind_fields: xr.Dataset = xr.open_dataset(wind_speeds_path)
    logging.info(wind_fields.max_wind_speed)  # use xarray repr

    logging.info(f"Using damage thresholds: {speed_thresholds} [m/s]")

    logging.info("Simulating electricity network failure due to wind damage...")

    args = ((storm_id, wind_fields.max_wind_speed, splits, speed_thresholds, network) for storm_id in wind_fields.event_id)

    exposure_by_storm: list[xr.Dataset] = []
    if parallel:
        with multiprocessing.Pool() as pool:
            exposure_by_storm = pool.starmap(degrade_grid_with_storm, args)
    else:
        for arg in args:
            exposure_by_storm.append(degrade_grid_with_storm(*arg))

    exposure = xr.merge(exposure_by_storm).to_netcdf(damages_path)
