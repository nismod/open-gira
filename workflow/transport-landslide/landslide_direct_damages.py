"""
Given an exposure estimate and some damage curves, calculate the damage
fraction for exposed assets.
"""

import logging
import sys
import warnings

import geopandas as gpd
import pandas as pd
from snail.damages import PiecewiseLinearDamageCurve

from open_gira import fields
from open_gira.direct_damages import annotate_rehab_cost
from open_gira.io import write_empty_frames, read_damage_curves, read_rehab_costs


if __name__ == "__main__":
    try:
        unsplit_path: str = snakemake.input["unsplit"]
        exposure_path: str = snakemake.input["exposure"]
        rehabilitation_costs_path: str = snakemake.input["rehab_cost"]
        damage_curves_dir: str = snakemake.input["damage_curves"]
        split_ead_and_cost_per_trigger_path: str = snakemake.output[
            "split_ead_and_cost_per_trigger"
        ]
        ead_and_cost_per_trigger_path: str = snakemake.output[
            "ead_and_cost_per_trigger"
        ]
        network_type: str = snakemake.params["network_type"].split("-")[0]
    except NameError:
        raise ValueError("Must be run via snakemake.")

    OUTPUT_FILE_PATHS: tuple[str] = (
        split_ead_and_cost_per_trigger_path,
        ead_and_cost_per_trigger_path,
    )
    HAZARD_TYPE = "landslide"

    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    # load curves first so if we fail here, we've failed early
    # and we don't try and load the (potentially large) exposure file
    damage_curves_all = read_damage_curves(
        damage_curves_dir, HAZARD_TYPE, set((network_type,))
    )
    damage_curve_data = damage_curves_all[network_type]
    assert (
        "occurrence" in damage_curve_data.columns
    ), "Expected 'occurrence' column in landslide damage curve"

    # Parse damage curve data into dict of DamageCurve objects
    damage_curves = {}
    for ratio_col in [c for c in damage_curve_data.columns if c != "occurrence"]:
        damage_curves[f"{network_type}_{ratio_col}"] = PiecewiseLinearDamageCurve(
            damage_curve_data[["occurrence", ratio_col]].rename(
                columns={"occurrence": "intensity", ratio_col: "damage"}
            )
        )
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
    rehab_cost = read_rehab_costs(rehabilitation_costs_path)
    exposure = annotate_rehab_cost(exposure, network_type, rehab_cost)

    # column groupings for data selection
    initial_hazard_columns = [
        col for col in exposure.columns if col.startswith(fields.HAZARD_PREFIX)
    ]
    exposure[f"{fields.HAZARD_PREFIX}_{HAZARD_TYPE}_sum"] = exposure[
        initial_hazard_columns
    ].sum(axis=1)

    hazard_columns = [
        col for col in exposure.columns if col.startswith(fields.HAZARD_PREFIX)
    ]
    non_hazard_columns = list(set(exposure.columns) - set(hazard_columns))

    #############################
    # EXPECTED ANNUAL DAMAGE COST
    #############################
    direct_damages = {}
    for hazard_probability_column in hazard_columns:
        # hazard maps give probability of occurrence
        hazard_probability = exposure[hazard_probability_column]
        # any non-zero probability of landslide has an "occurrence" value of 1.0
        exposure_intensity = (hazard_probability > 0).astype("float")
        for damage_curve_key, damage_curve in damage_curves.items():
            # damage curves are step functions based on 0-1 occurrence
            damage_fraction = damage_curve.damage_fraction(exposure_intensity)
            # damage cost is calculated directly from damage fraction
            damage_cost = damage_fraction * exposure[fields.REHAB_COST]
            # and so expected damage is (exposed value * damage fraction * probability of occurrence)
            expected_damage = damage_cost * hazard_probability
            direct_damages[f"{hazard_probability_column}__{damage_curve_key}_EAD"] = (
                expected_damage
            )
            direct_damages[
                f"{hazard_probability_column}__{damage_curve_key}_damage_cost"
            ] = damage_cost

    direct_damages = pd.DataFrame(direct_damages)
    split_ead_and_cost_per_trigger = pd.concat(
        [exposure[non_hazard_columns], direct_damages], axis=1
    )
    grouped_direct_damages = (
        pd.concat([exposure[["id"]], direct_damages], axis=1).groupby("id").sum()
    )

    #########################################
    # JOINING, VALIDATION AND SERIALIZATION #
    #########################################

    logging.info("Reading raw network data for unsplit geometry")
    unsplit: gpd.GeoDataFrame = gpd.read_parquet(unsplit_path)
    logging.info(f"{unsplit.shape=}")

    # lose columns like "cell_indices" or rastered length measures that are specific to _rastered_ edges
    non_hazard_output_columns = list(set(non_hazard_columns) & set(unsplit.columns))
    unsplit_subset = unsplit[non_hazard_output_columns].set_index("id", drop=False)

    # rejoin direct damage cost estimates with geometry and metadata columns and write to disk
    # join on 'right' / grouped_direct_damages index to only keep rows we have damages for
    ead_and_cost_per_trigger = unsplit_subset.join(
        grouped_direct_damages, validate="one_to_one", how="right"
    )
    # we may not have calculated damages for every possible asset_type
    assert len(ead_and_cost_per_trigger) <= len(unsplit_subset)
    assert "id" in ead_and_cost_per_trigger.columns

    logging.info(
        f"Writing out {split_ead_and_cost_per_trigger.shape=} "
        "(per unified geometry, hazard RP map and hazard map (integrated RP))"
    )
    split_ead_and_cost_per_trigger.to_parquet(split_ead_and_cost_per_trigger_path)

    logging.info(
        f"Writing out {ead_and_cost_per_trigger.shape=} "
        "(per unified geometry, hazard RP map and hazard map (integrated RP))"
    )
    ead_and_cost_per_trigger.to_parquet(ead_and_cost_per_trigger_path)

    logging.info("Done calculating direct damages")
