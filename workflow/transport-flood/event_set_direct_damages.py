"""
Given an exposure estimate and some damage curves, calculate the damage
fraction for exposed assets.
"""

import logging
import sys
import warnings

import geopandas as gpd

from open_gira import fields
from open_gira.direct_damages import annotate_rehab_cost, direct_damage
from open_gira.io import write_empty_frames, read_damage_curves, read_rehab_costs


if __name__ == "__main__":

    try:
        unsplit_path: str = snakemake.input["unsplit"]
        exposure_path: str = snakemake.input["exposure"]
        rehabilitation_costs_path: str = snakemake.input["rehab_cost"]
        damage_curves_dir: str = snakemake.input["damage_curves"]
        damage_fraction_path: str = snakemake.output["damage_fraction"]
        damage_cost_path: str = snakemake.output["damage_cost"]
        network_type: str = snakemake.params["network_type"]
        hazard_type: str = snakemake.params["hazard_type"]
        # note, this config entry might have been mutated on execution of the Snakefile
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
    rehab_cost = read_rehab_costs(rehabilitation_costs_path)
    exposure = annotate_rehab_cost(exposure, network_type, rehab_cost)

    # column groupings for data selection
    hazard_columns = [col for col in exposure.columns if col.startswith(fields.HAZARD_PREFIX)]
    non_hazard_columns = list(set(exposure.columns) - set(hazard_columns))

    damage_fraction, grouped_direct_damages = direct_damage(exposure, damage_curves, hazard_columns, non_hazard_columns)

    logging.info("Reading raw network data for unsplit geometry")
    unsplit: gpd.GeoDataFrame = gpd.read_parquet(unsplit_path)
    logging.info(f"{unsplit.shape=}")

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