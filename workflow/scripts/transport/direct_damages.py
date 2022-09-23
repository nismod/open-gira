"""
Given an exposure estimate and some damage curves, calculate the damage ratio
for exposed assets.
"""

import logging
import sys
import warnings

import pandas as pd

import utils


if __name__ == "__main__":
    try:
        exposure_path: str = snakemake.input["exposure"]
        damage_fraction_path: str = snakemake.output["damage_fraction"]
        damage_curves_path: str = snakemake.config["direct_damages"]["curves_path"]
        network_type: str = snakemake.params["network_type"]
        hazard_type: str = snakemake.params["hazard_type"]
        asset_types: set[str] = set(snakemake.config["direct_damages"]["asset_types"])
    except NameError:
        # running as a script, rather than via snakemake...
        (
            exposure_path,
            damage_fraction_path,
            damage_curves_path,
            network_type,
            hazard_type,
            asset_types
        ) = sys.argv[1:]
        # where for example:
        # exposure_path = results/splits/egypt-latest_filter-road/hazard-aqueduct-river/slice-1.parquet
        # damage_fraction_path = results/direct_damages/egypt-latest_filter-road/hazard-aqueduct-river/slice-1.parquet
        # damage_curves_path = bundled_data/damage_curves.xlsx
        # network_type = road
        # hazard_type = flood
        # asset_types = "road_unpaved,road_paved,road_bridge"

        # parse comma seperated string into set of strings
        asset_types: set = {s.strip() for s in asset_types.split(",")}

    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    exposure: pd.DataFrame = pd.read_parquet(exposure_path)

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
    hazard_columns = [col for col in exposure.columns if col.startswith("hazard-")]
    non_hazard_columns = list(set(exposure.columns) - set(hazard_columns))

    # calculate damages
    direct_damages_by_asset_type = []
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

        asset_type_mask = exposure.asset_type == asset_type
        asset_exposure = exposure.loc[asset_type_mask, hazard_columns]

        # apply damage_curve function element-wise
        # TODO: Vectorise this (take a series of intensities and return a
        # series of damage fractions.
        direct_damages = asset_exposure.applymap(damage_curve)

        # store the computed direct damages and any columns we started with
        # other than exposure
        direct_damages_by_asset_type.append(
            direct_damages.join(exposure.loc[asset_type_mask, non_hazard_columns])
        )

    # concatenate into single dataframe and save to disk
    pd.concat(direct_damages_by_asset_type).to_parquet(damage_fraction_path)
