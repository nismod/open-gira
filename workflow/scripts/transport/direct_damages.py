"""
Given an exposure estimate and some damage curves, calculate the damage
fraction for exposed assets.
"""

import logging
import sys
import warnings
from glob import glob
from os.path import splitext, basename, relpath, join

import geopandas as gpd
import pandas as pd
from scipy.interpolate import interp1d

import utils
from plot_damage_distributions import natural_sort


# exposure table hazard intensity fields expected to be prefixed as such
HAZARD_INTENSITY_PREFIX = "hazard-"
# exposure table field containing the cost to rebuild per unit length
REHABILITATION_COST_FIELD = "rehab_cost_USD_per_km"
# length of edge calculated after intersection
SPLIT_LENGTH_FIELD = "length_km"


def load_damage_curves(damage_curves_dir: str, hazard_type: str, asset_types: set) -> dict[str, pd.DataFrame]:
    """
    Load damage curves from CSVs on disk

    Expected to reside in following structure:
    <damage_curves_dir>/<hazard_type>/<asset_type>.csv

    Damage curve files may have comments, these are lines starting with COMMENT_PREFIX

    Args:
        damage_curves_dir (str): Path to folder containing hazards
        hazard_type (str): Name of hazard folder containing asset specific curves
        asset_types (set): Asset types we require damage curves for

    Returns (dict[str, pd.DataFrame):
        Mapping from asset_type to respective damage curve
    """

    # lines beginning with this character will be ignored by pandas
    COMMENT_PREFIX: str = "#"

    # fetch damage curves for relevant hazard type
    damage_curve_paths = glob(join(damage_curves_dir, hazard_type, "*.csv"))

    damage_curves: dict[str, pd.DataFrame] = {
        # curves expected to be named as a value of Asset class, e.g. RoadAssets.BRIDGE -> road_bridge.csv
        # dict is asset_type: dataframe with hazard intensity [0, inf] and damage fraction [0, 1]
        splitext(basename(path))[0]: pd.read_csv(path, comment=COMMENT_PREFIX) for path in damage_curve_paths
    }

    for asset_type, damage_curve in damage_curves.items():
        # check hazard intensity and damage fraction values are 0 or positive real
        assert ((damage_curve >= 0).all()).all()
        # check damage fraction is less than or equal to 1
        assert (damage_curve.iloc[:, 1] <= 1).all()

    if not set(damage_curves.keys()).issuperset(asset_types):
        raise RuntimeError(f"requested {asset_types=} not all found: {damage_curves.keys()=}")

    return damage_curves


if __name__ == "__main__":

    try:
        unsplit_path: str = snakemake.input["unsplit"]
        exposure_path: str = snakemake.input["exposure"]
        damage_fraction_path: str = snakemake.output["damage_fraction"]
        damage_cost_path: str = snakemake.output["damage_cost"]
        damage_curves_dir: str = snakemake.config["direct_damages"]["curves_dir"]
        network_type: str = snakemake.params["network_type"]
        hazard_type: str = snakemake.params["hazard_type"]
        asset_types: set[str] = set(snakemake.config["direct_damages"]["asset_types"])

    except NameError:
        # running as a script, rather than via snakemake...
        (
            unsplit_path,
            exposure_path,
            damage_fraction_path,
            damage_cost_path,
            damage_curves_dir,
            network_type,
            hazard_type,
            asset_types
        ) = sys.argv[1:]

        # where for example:
        # unsplit_path = results/geoparquet/egypt-latest_filter-road/slice-0.geoparquet
        # exposure_path = results/splits/egypt-latest_filter-road/hazard-aqueduct-river/slice-1.geoparquet
        # damage_fraction_path = results/direct_damages/egypt-latest_filter-road/hazard-aqueduct-river/slice-1_fraction.parquet
        # damage_cost_path = results/direct_damages/egypt-latest_filter-road/hazard-aqueduct-river/slice-1_cost.parquet
        # damage_curves_dir = bundled_data/damage_curves/
        # network_type = road
        # hazard_type = flood
        # asset_types = "road_unpaved,road_paved,road_bridge"

        # parse comma separated string into set of strings
        asset_types: set = {s.strip() for s in asset_types.split(",")}

    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    # load curves first so if we fail here, we've failed early
    # and we don't try and load the (potentially large) exposure file
    damage_curves = load_damage_curves(damage_curves_dir, hazard_type, asset_types)
    logging.info(f"Available damage curves: {damage_curves.keys()}")

    logging.info(f"Reading exposure (network/raster intersection) data")
    exposure: gpd.GeoDataFrame = gpd.read_parquet(exposure_path)
    logging.info(f"{exposure.shape=}")

    if exposure.empty:
        logging.info("No data in geometry column, writing empty files.")

        # snakemake requires that output files exist, even if empty
        for path in (damage_fraction_path, damage_cost_path):
            utils.write_empty_frames(path)
        sys.exit(0)  # exit gracefully so snakemake will continue

    # column groupings for data selection
    hazard_columns = [col for col in exposure.columns if col.startswith(HAZARD_INTENSITY_PREFIX)]
    non_hazard_columns = list(set(exposure.columns) - set(hazard_columns))

    # calculate damages for assets we have damage curves for
    damage_fraction_by_asset_type = []
    logging.info(f"Exposed assets {set(exposure.asset_type)}")
    for asset_type in set(exposure.asset_type) & set(damage_curves.keys()):

        logging.info(f"Processing {asset_type=}")
        damage_curve: pd.DataFrame = damage_curves[asset_type]

        # pick out rows of asset type and columns of hazard intensity
        asset_type_mask: gpd.GeoDataFrame = exposure.asset_type == asset_type
        asset_exposure: pd.DataFrame = pd.DataFrame(exposure.loc[asset_type_mask, hazard_columns])

        # create interpolated damage curve for given asset type
        hazard_intensity, damage_fraction = damage_curve.iloc[:, 0], damage_curve.iloc[:, 1]
        # if curve of length n, where x < x_0, y = y_0 and where x > x_n, y = y_n
        bounds = tuple(f(damage_fraction) for f in (min, max))
        interpolated_damage_curve = interp1d(hazard_intensity, damage_fraction, kind='linear', fill_value=bounds, bounds_error=False)

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
    # write to disk
    logging.info(f"Writing {damage_fraction.shape=} to disk")
    damage_fraction.to_parquet(damage_fraction_path)

    # multiply the damage fraction estimates by a cost to rebuild the asset
    # units are: 1 * USD/km * km = USD
    logging.info("Calculating direct damage costs")
    direct_damages_only = damage_fraction[hazard_columns] \
            .multiply(damage_fraction[REHABILITATION_COST_FIELD], axis="index") \
            .multiply(damage_fraction[SPLIT_LENGTH_FIELD], axis="index")

    logging.info(f"Reading raw network data for unsplit geometry")
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

    # lose columns like "cell_indicies" or rastered length measures that are specific to _rastered_ edges
    non_hazard_output_columns = list(set(non_hazard_columns) & set(unsplit.columns))
    unsplit_subset = unsplit[non_hazard_output_columns].set_index("edge_id")

    # join the direct damage estimates with the unsplit geometry and other relevant columns
    direct_damages = unsplit_subset.join(grouped_direct_damages, validate="one_to_one")
    direct_damages["edge_id"] = direct_damages.index

    # sort damage columns alphabetically
    direct_damages = direct_damages[natural_sort(direct_damages.columns)]

    assert len(unsplit) == len(direct_damages)

    # write to disk
    logging.info(f"Writing {direct_damages.shape=} to disk")
    direct_damages.to_parquet(damage_cost_path)

    logging.info("Done")
