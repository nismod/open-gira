"""
For each maximum wind field associated with a storm:
    1) Fail edge segments who experience a wind speed greater than given threshold
    2) Attempt to allocate power from sources to sinks over degraded network
    3) Calculate ratio of nominal power to degraded power, the 'supply factor'
"""

import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

from open_gira import fields
from open_gira.grid import allocate_power_to_targets
import snkit


logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)


def main():
    edges_path: str = snakemake.input.grid_edges
    nodes_path: str = snakemake.input.grid_nodes
    splits_path: str = snakemake.input.grid_splits
    wind_speeds_path: str = snakemake.input.wind_speeds
    # sort into ascending order; if no damage at a given threshold,
    # more resilient thresholds are guaranteed to be safe
    speed_thresholds: list[float] = sorted(snakemake.config["transmission_windspeed_failure"])
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

    speeds: xr.DataSet = xr.open_dataset(wind_speeds_path)
    logging.info(speeds)  # use xarray repr

    logging.info(f"Using damage thresholds: {speed_thresholds} [m/s]")

    logging.info("Simulating electricity network failure due to wind damage...")

    # build object for accumulating grid damage results
    target_ids = network.nodes[network.nodes.asset_type == "target"].id.values
    coords_shape = (len(speeds.event_id), len(speed_thresholds), len(target_ids))
    empty_cube = np.full(coords_shape, np.nan)
    dim_names = ["event_id", "threshold", "target"]
    ds = xr.Dataset(
        data_vars=dict(
            supply_factor=(dim_names, empty_cube),
            customers_affected=(dim_names, empty_cube)
        ),
        coords=dict(
            event_id=speeds.event_id,
            threshold=speed_thresholds,
            target=target_ids
        )
    )

    for storm_id in speeds.event_id:

        # rank 1, length of splits DataFrame
        # N.B. to index at points rather than the cross-product of indicies, index with DataArrays
        # https://docs.xarray.dev/en/stable/user-guide/indexing.html#vectorized-indexing
        max_wind_speeds: xr.DataArray = speeds.max_wind_speed.sel(event_id=storm_id).isel(
            long=splits[fields.RASTER_I].to_xarray(),
            lat=splits[fields.RASTER_J].to_xarray()
        )

        storm_id = str(storm_id.values)  # coerce singleton array into str
        logging.info(storm_id)

        for threshold in speed_thresholds:
            survival_mask: pd.Series = (max_wind_speeds < threshold).to_pandas()

            try:
                n_failed: int = survival_mask.value_counts()[False]
            except KeyError:
                break  # network intact

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
            failure_str = "-> {:>6.2f}% edges failed @ {:.1f} [m/s] threshold, {:d} -> {:d} components"
            logging.info(failure_str.format(100 * fraction_failed, threshold, c_nominal, c_surviving))

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
            ds.supply_factor.loc[indicies] = targets.supply_factor
            ds.customers_affected.loc[indicies] = targets.customers_affected

    ds.to_netcdf(damages_path)


if __name__ == "__main__":
    main()
