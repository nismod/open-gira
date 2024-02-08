"""
Given an exposure estimate and some damage curves, calculate the damage
fraction for exposed assets.
"""

import logging
import sys
import warnings

import geopandas as gpd
import pandas as pd
from scipy.interpolate import interp1d

from open_gira import fields
from open_gira.direct_damages import annotate_rehab_cost
from open_gira.io import write_empty_frames, read_damage_curves
from open_gira.utils import natural_sort


if __name__ == "__main__":

    try:
        unsplit_path: str = snakemake.input["unsplit"]
        exposure_path: str = snakemake.input["exposure"]
        damage_fraction_path: str = snakemake.output["damage_fraction"]
        damage_cost_path: str = snakemake.output["damage_cost"]
        damage_curves_dir: str = snakemake.config["direct_damages"]["curves_dir"]
        rehabilitation_costs_path = snakemake.config["transport"]["rehabilitation_costs_path"]
        network_type: str = snakemake.params["network_type"]
        hazard_type: str = snakemake.params["hazard_type"]
        asset_types: set[str] = set(snakemake.config["direct_damages"]["asset_types"])
    except NameError:
        raise ValueError("Must be run via snakemake.")

    OUTPUT_FILE_PATHS: tuple[str] = (
        damage_fraction_path,
        damage_cost_path,
    )

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    # load curves first so if we fail here, we've failed early
    # and we don't try and load the (potentially large) exposure file
    damage_curves = read_damage_curves(damage_curves_dir, hazard_type, asset_types)
    logging.info(f"Available damage curves: {damage_curves.keys()}")

    logging.info("Reading exposure (network/raster intersection) data")
    exposure: gpd.GeoDataFrame = gpd.read_parquet(exposure_path)
    logging.info(f"{exposure.shape=}")

    if exposure.empty:
        logging.info("No data in geometry column, writing empty files.")

        # snakemake requires that output files exist, even if empty
        for path in OUTPUT_FILE_PATHS:
            write_empty_frames(path)
        sys.exit(0)  # exit gracefully so snakemake will continue

    logging.info("Annotate network with rehabilitation costs")
    rehab_cost = pd.read_excel(rehabilitation_costs_path, sheet_name=network_type)
    exposure = annotate_rehab_cost(exposure, network_type, rehab_cost)

    # column groupings for data selection
    hazard_columns = [col for col in exposure.columns if col.startswith(fields.HAZARD_PREFIX)]
    non_hazard_columns = list(set(exposure.columns) - set(hazard_columns))

    ##########################################################
    # DAMAGE FRACTIONS (per split geometry, for all rasters) #
    ##########################################################

    # calculate damages for assets we have damage curves for
    damage_fraction_by_asset_type = []
    logging.info(f"Exposed assets {set(exposure.asset_type)}")
    for asset_type in natural_sort(set(exposure.asset_type) & set(damage_curves.keys())):

        logging.info(f"Processing {asset_type=}")
        damage_curve: pd.DataFrame = damage_curves[asset_type]

        # pick out rows of asset type and columns of hazard intensity
        asset_type_mask: gpd.GeoDataFrame = exposure.asset_type == asset_type
        asset_exposure: pd.DataFrame = pd.DataFrame(exposure.loc[asset_type_mask, hazard_columns])

        # create interpolated damage curve for given asset type
        hazard_intensity, damage_fraction = damage_curve.iloc[:, 0], damage_curve.iloc[:, 1]
        # if curve of length n, where x < x_0, y = y_0 and where x > x_n, y = y_n
        bounds = tuple(f(damage_fraction) for f in (min, max))
        interpolated_damage_curve = interp1d(
            hazard_intensity,
            damage_fraction,
            kind='linear',
            fill_value=bounds,
            bounds_error=False
        )

        # apply damage_curve function to exposure table
        # the return value of interpolated_damage_curve is a numpy array
        logging.info("Calculating damage fractions")
        damage_fraction_for_asset_type = pd.DataFrame(
            interpolated_damage_curve(asset_exposure),
            index=asset_exposure.index,
            columns=asset_exposure.columns
        )

        # store the computed direct damages and any columns we started with
        # (other than exposure)
        damage_fraction_by_asset_type.append(
            pd.concat(
                [
                    damage_fraction_for_asset_type,
                    exposure.loc[asset_type_mask, non_hazard_columns]
                ],
                axis="columns"
            )
        )

    # concatenate damage fractions for different asset types into single dataframe
    damage_fraction: gpd.GeoDataFrame = gpd.GeoDataFrame(pd.concat(damage_fraction_by_asset_type))

    ###################################################################
    # DAMAGE COST (for split, then grouped geometry, for all rasters) #
    ###################################################################

    # multiply the damage fraction estimates by a cost to rebuild the asset
    # units are: 1 * USD/km * km = USD
    logging.info("Calculating direct damage costs")
    direct_damages_only = damage_fraction[hazard_columns] \
        .multiply(damage_fraction[fields.REHAB_COST], axis="index") \
        .multiply(damage_fraction[fields.SPLIT_LENGTH], axis="index")

    logging.info("Reading raw network data for unsplit geometry")
    unsplit: gpd.GeoDataFrame = gpd.read_parquet(unsplit_path)
    logging.info(f"{unsplit.shape=}")

    # join the other fields with the direct damage estimates
    logging.info("Unifying rasterised segments and summing damage costs")

    # grouping on edge_id, sum all direct damage estimates to give a total dollar cost per edge
    direct_damages = pd.concat(
        [direct_damages_only, damage_fraction["edge_id"]],
        axis="columns"
    ).set_index("edge_id")
    grouped_direct_damages = direct_damages.groupby(direct_damages.index).sum()

    #########################################
    # JOINING, VALIDATION AND SERIALIZATION #
    #########################################

    # lose columns like "cell_indicies" or rastered length measures that are specific to _rastered_ edges
    non_hazard_output_columns = list(set(non_hazard_columns) & set(unsplit.columns))
    unsplit_subset = unsplit[non_hazard_output_columns].set_index("edge_id", drop=False)

    # rejoin direct damage cost estimates with geometry and metadata columns and write to disk
    # join on 'right' / grouped_direct_damages index to only keep rows we have damages for
    direct_damages = unsplit_subset.join(grouped_direct_damages, validate="one_to_one", how="right")
    direct_damages["edge_id"] = direct_damages.index
    # we may not have calculated damages for every possible asset_type
    assert len(direct_damages) <= len(unsplit_subset)
    assert "edge_id" in direct_damages.columns

    # damage_fraction is on the split geometries, will have more rows
    assert len(damage_fraction) >= len(direct_damages)

    for dataframe in (damage_fraction, direct_damages):
        assert "edge_id" in dataframe

    logging.info(f"Writing out {damage_fraction.shape=} (per split geometry, event raster)")
    damage_fraction.to_parquet(damage_fraction_path)

    logging.info(f"Writing out {direct_damages.shape=} (per unified geometry, event raster)")
    direct_damages.to_parquet(damage_cost_path)

    logging.info("Done calculating direct damages")
