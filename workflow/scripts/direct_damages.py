"""
Given an exposure estimate and some damage curves, calculate the damage
fraction for exposed assets.
"""

from collections import defaultdict
import logging
import sys
import warnings

import geopandas as gpd
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import simpson

from open_gira import fields
from open_gira.direct_damages import ReturnPeriodMap, generate_rp_maps
from open_gira.io import write_empty_frames, read_damage_curves
from open_gira.utils import natural_sort


if __name__ == "__main__":

    try:
        unsplit_path: str = snakemake.input["unsplit"]
        exposure_path: str = snakemake.input["exposure"]
        damage_fraction_path: str = snakemake.output["damage_fraction"]
        damage_cost_path: str = snakemake.output["damage_cost"]
        expected_annual_damages_path: str = snakemake.output["expected_annual_damages"]
        return_period_and_ead_path: str = snakemake.output["return_period_and_ead"]
        damage_curves_dir: str = snakemake.config["direct_damages"]["curves_dir"]
        network_type: str = snakemake.params["network_type"]
        hazard_type: str = snakemake.params["hazard_type"]
        asset_types: set[str] = set(snakemake.config["direct_damages"]["asset_types"])
    except NameError:
        raise ValueError("Must be run via snakemake.")

    OUTPUT_FILE_PATHS: tuple[str] = (
        damage_fraction_path,
        damage_cost_path,
        expected_annual_damages_path,
        return_period_and_ead_path
    )

    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

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

    ###############################################################################
    # EXPECTED ANNUAL DAMAGE COST (for grouped geometry, aggregations of rasters) #
    ###############################################################################

    # reduce climate model / subsidence to simple MIN/MAX
    # generate a mapping from a 'family' of hazards to their set of related return period maps
    model_families: dict[str, set[ReturnPeriodMap]] = defaultdict(set)
    for rp_map in generate_rp_maps(grouped_direct_damages.columns, prefix=fields.HAZARD_PREFIX):
        model_families[rp_map.without_model].add(rp_map)  # only differ by climate model / subsidence

    aggregations = ("min", "mean", "max")
    logging.info(f"Applying {aggregations=} to input raster models")
    for family_name, family_rp_maps in model_families.items():
        for agg_str in aggregations:
            sample_map, *_ = family_rp_maps
            family_aggregation_name = sample_map.name.replace(sample_map.model, agg_str.upper())
            family_column_names: list[str] = [f"{fields.HAZARD_PREFIX}{rp_map.name}" for rp_map in family_rp_maps]
            agg_func = getattr(grouped_direct_damages[family_column_names], agg_str)
            grouped_direct_damages[f"{fields.HAZARD_PREFIX}{family_aggregation_name}"] = agg_func(axis="columns")
        grouped_direct_damages = grouped_direct_damages.drop(columns=family_column_names)

    # integrate over return periods for expected annual damages
    rp_map_families: dict[str, set[ReturnPeriodMap]] = defaultdict(set)
    for rp_map in generate_rp_maps(grouped_direct_damages.columns, prefix=fields.HAZARD_PREFIX):
        rp_map_families[rp_map.without_RP].add(rp_map)  # only differ by return period

    expected_annual_damages = {}
    logging.info(f"Integrating {len(rp_map_families)} damage-probability curves")
    for family_name, family_rp_maps in rp_map_families.items():

        # sort by least to most probable
        sorted_rp_maps: list[ReturnPeriodMap] = sorted(family_rp_maps)

        # [0, 1] valued decimal probabilities
        probabilities: list[float] = [rp_map.annual_probability for rp_map in sorted_rp_maps]
        # family subset of grouped_direct_damages
        family_column_names: list[str] = [f"{fields.HAZARD_PREFIX}{rp_map.name}" for rp_map in sorted_rp_maps]
        family_direct_damages: pd.DataFrame = grouped_direct_damages[family_column_names]

        # integrate the damage as a function of probability curve using Simpson's rule
        # Simpson's rule as the function to be integrated is non-linear
        # add _EAD for easier downstream selection of these columns vs. RPs
        expected_annual_damages[f"{fields.HAZARD_PREFIX}{family_name}_EAD"] = \
            simpson(family_direct_damages, x=probabilities, axis=1)

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

    expected_annual_damages_only = pd.DataFrame(data=expected_annual_damages, index=grouped_direct_damages.index)
    # rejoin expected annual damage cost estimates with geometry and metadata columns and write to disk
    # join on 'right' / expected_annual_damages index to only keep rows we have damages for
    expected_annual_damages = gpd.GeoDataFrame(
        unsplit_subset.join(expected_annual_damages_only, validate="one_to_one", how="right")
    )
    assert len(expected_annual_damages) <= len(unsplit_subset)
    assert "edge_id" in expected_annual_damages.columns

    # combined the per return period and the integrated outputs into a single dataframe
    return_period_and_ead_damages = direct_damages.join(expected_annual_damages_only, validate="one_to_one")
    assert len(return_period_and_ead_damages) == len(direct_damages) == len(expected_annual_damages_only)

    # damage_fraction is on the split geometries, will have more rows
    assert len(damage_fraction) >= len(direct_damages)
    assert len(damage_fraction) >= len(expected_annual_damages)

    # direct_damages and expected_annual_damages should have the same feature count
    assert len(direct_damages) == len(expected_annual_damages)

    for dataframe in (damage_fraction, direct_damages, expected_annual_damages, return_period_and_ead_damages):
        assert "edge_id" in dataframe

    logging.info(f"Writing out {damage_fraction.shape=} (per split geometry, hazard RP map)")
    damage_fraction.to_parquet(damage_fraction_path)

    logging.info(f"Writing out {direct_damages.shape=} (per unified geometry, hazard RP map)")
    direct_damages.to_parquet(damage_cost_path)

    logging.info(f"Writing out {expected_annual_damages.shape=} (per unified geometry, hazard map (integrated RP))")
    expected_annual_damages.to_parquet(expected_annual_damages_path)

    logging.info(
        f"Writing out {return_period_and_ead_damages.shape=} "
        "(per split geometry, hazard RP map and hazard map (integrated RP))"
    )
    return_period_and_ead_damages.to_parquet(return_period_and_ead_path)

    logging.info("Done calculating direct damages")