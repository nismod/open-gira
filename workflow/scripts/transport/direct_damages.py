"""
Given an exposure estimate and some damage curves, calculate the damage
fraction for exposed assets.
"""

import logging
import sys
import warnings

import geopandas as gpd
import pandas as pd

import utils


# exposure table hazard intensity fields expected to be prefixed as such
HAZARD_INTENSITY_PREFIX = "hazard-"
# exposure table field containing the cost to rebuild per unit length
REHABILITATION_COST_FIELD = "rehab_cost_USD_per_km"


if __name__ == "__main__":

    try:
        exposure_path: str = snakemake.input["exposure"]
        damage_fraction_path: str = snakemake.output["damage_fraction"]
        damage_cost_path: str = snakemake.output["damage_cost"]
        damage_curves_path: str = snakemake.config["direct_damages"]["curves_path"]
        network_type: str = snakemake.params["network_type"]
        hazard_type: str = snakemake.params["hazard_type"]
        asset_types: set[str] = set(snakemake.config["direct_damages"]["asset_types"])

    except NameError:
        # running as a script, rather than via snakemake...
        (
            exposure_path,
            damage_fraction_path,
            damage_cost_path,
            damage_curves_path,
            network_type,
            hazard_type,
            asset_types
        ) = sys.argv[1:]

        # where for example:
        # exposure_path = results/splits/egypt-latest_filter-road/hazard-aqueduct-river/slice-1.parquet
        # damage_fraction_path = results/direct_damages/egypt-latest_filter-road/hazard-aqueduct-river/slice-1_fraction.parquet
        # damage_cost_path = results/direct_damages/egypt-latest_filter-road/hazard-aqueduct-river/slice-1_cost.parquet
        # damage_curves_path = bundled_data/damage_curves.xlsx
        # network_type = road
        # hazard_type = flood
        # asset_types = "road_unpaved,road_paved,road_bridge"

        # parse comma separated string into set of strings
        asset_types: set = {s.strip() for s in asset_types.split(",")}

    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    exposure: gpd.GeoDataFrame = gpd.read_parquet(exposure_path)

    if exposure.empty:
        logging.info("writing empty files and skipping processing...")

        # snakemake requires that output files exist, even if empty
        utils.write_empty_frames(damage_fraction_path)
        sys.exit(0)  # exit gracefully so snakemake will continue

    # fetch damage curves for relevant assets
    damage_curves: pd.DataFrame = pd.read_excel(damage_curves_path, sheet_name=hazard_type)
    asset_types_with_damage_curves = asset_types & set(damage_curves.columns)
    if len(asset_types_with_damage_curves) == 0:
        raise RuntimeError(f"none of {asset_types=} found in {damage_curves.columns=}")

    logging.info(f"Found {asset_types_with_damage_curves=}")

    # column groupings for data selection
    hazard_columns = [col for col in exposure.columns if col.startswith(HAZARD_INTENSITY_PREFIX)]
    non_hazard_columns = list(set(exposure.columns) - set(hazard_columns))

    # calculate damages
    damage_fraction_by_asset_type = []
    for asset_type in asset_types_with_damage_curves:

        def damage_curve(i: float) -> float:
            """
            Inner function acting as a damage curve for a given asset type.

            Given a hazard intensity scalar, return the damage fraction. If the
            intensity is outside the range provided in `damage_curves`, clip to
            the minimum or maximum value as appropriate.

            Arguments:
                i (float): Hazard intensity

            Returns:
                float: Damage fraction. N.B. Expected to be in range [0, 1].
            """

            # we assume the first col in the table is the hazard intensity array
            nearest_intensity_idx = (damage_curves.iloc[:, 0] - i).abs().idxmin()

            return damage_curves.loc[nearest_intensity_idx, asset_type]

        asset_type_mask: gpd.GeoDataFrame = exposure.asset_type == asset_type
        asset_exposure: pd.DataFrame = pd.DataFrame(exposure.loc[asset_type_mask, hazard_columns])

        # apply damage_curve function element-wise
        # TODO: vectorise this?
        damage_fraction_for_asset_type: pd.DataFrame = asset_exposure.applymap(damage_curve)

        # store the computed direct damages and any columns we started with
        # (other than exposure)
        damage_fraction_by_asset_type.append(
            damage_fraction_for_asset_type.join(
                exposure.loc[asset_type_mask, non_hazard_columns]
            )
        )

    # concatenate damage fractions for different asset types into single dataframe
    damage_fraction: gpd.GeoDataFrame = gpd.GeoDataFrame(pd.concat(damage_fraction_by_asset_type))
    # write to disk
    damage_fraction.to_parquet(damage_fraction_path)

    # multiply the damage fraction estimates by a cost to rebuild the asset
    # units are: 1 * USD/km * km = USD
    direct_damages_only: pd.DataFrame = pd.DataFrame(
        damage_fraction.loc[:, hazard_columns] \
            .multiply(exposure.loc[:, REHABILITATION_COST_FIELD], axis="index") \
            .multiply(exposure.loc[:, "length_km"], axis="index")
    )

    # join the other fields with the direct damage estimates
    direct_damages: gpd.GeoDataFrame = gpd.GeoDataFrame(
        direct_damages_only.join(exposure.loc[:, non_hazard_columns])
    )
    # write to disk
    direct_damages.to_parquet(damage_cost_path)
