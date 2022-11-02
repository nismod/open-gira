"""
Given an exposure estimate and some damage curves, calculate the damage
fraction for exposed assets.
"""

from abc import ABC, abstractmethod
from collections import defaultdict
from glob import glob
import logging
from os.path import splitext, basename, join
import re
import sys
from typing import Union
import warnings

import geopandas as gpd
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import simpson

import utils
from plot_damage_distributions import natural_sort


# exposure table hazard intensity fields expected to be prefixed as such
HAZARD_PREFIX = "hazard-"
# exposure table field containing the cost to rebuild per unit length
REHABILITATION_COST_FIELD = "rehab_cost_USD_per_km"
# length of edge calculated after intersection
SPLIT_LENGTH_FIELD = "length_km"


class ReturnPeriodMap(ABC):
    """
    Abstract class defining interface for return period hazard maps.
    """

    # identifying string containing other, inferred attributes, should be
    # unique among any collection of maps
    name: str
    # name of scenario, e.g. rcp4p5, historical
    scenario: str
    # year for which hazard map is valid (may be past, present or future)
    year: int
    # expect hazard to recur on average every return_period years
    return_period_years: float

    def __init__(self):
        if type(self) is ReturnPeriodMap:
            raise RuntimeError("Abstract class; please subclass rather than directly instantiating")

    @property
    @abstractmethod
    def without_RP(self):
        """
        A name identifying the kind of hazard, without any return period.

        N.B. For collapsing return periods into expected annual damages (EAD)
        it is useful to generate a name without return period information.
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def without_model(self) -> str:
        """
        A name identifying the attributes of the hazard, without any model
        information (climate model / subsidence).
        """
        raise NotImplementedError

    @property
    def annual_probability(self) -> float:
        """Likelihood of hazard occurring in any given year"""
        return 1 / self.return_period_years

    def __eq__(self, other):
        """Maps are equal if their name is"""
        return self.name == other.name

    def __lt__(self, other):
        """
        Permit sorting by annual probability (i.e. least likely events first),
        the usual order for an expected annual damages integration.
        """
        return self.annual_probability < other.annual_probability

    def __hash__(self):
        """Map name string should capture all that is unique about object."""
        return hash(self.name)


class AqueductFlood(ReturnPeriodMap):
    """
    Class holding information about aqueduct return period flood maps.

    Each point in these raster flood maps is an inundation depth for a given
    combination of e.g. scenaro, climate model, year, return period
    (probability).
    """

    # there are two subcategories of aqueduct flood map
    COASTAL = "inuncoast"
    RIVERINE = "inunriver"

    # coastal models may or may not include a subsidence component
    WITH_SUBSIDENCE = "wtsub"
    WITHOUT_SUBSIDENCE = "nosub"

    def __init__(self, name: str):
        """
        Infer attributes from name.

        Arguments:
            name (str): Name string expected to be in one of the following formats:
                Riverine:
                    inunriver_rcp8p5_00000NorESM1-M_2080_rp00005
                Coastal:
                    inuncoast_rcp8p5_wtsub_2080_rp1000_0_perc_50
        """

        if len(name.split(".")) > 1:
            raise ValueError(f"{name=} contains dots; remove any file extension")

        # store the original string for later use
        self.name = name

        map_type, *split_name = name.split("_")

        if map_type == self.RIVERINE:

            # unpack rest of name
            scenario, climate_model, year, return_period_years = split_name

            self.riverine = True
            self.coastal = False
            self.model = climate_model

        elif map_type == self.COASTAL:

            # unpack rest of name
            scenario, sub_str, year, return_period_years, *slr_perc_list = split_name

            if sub_str == self.WITH_SUBSIDENCE:
                subsidence = True
            elif sub_str == self.WITHOUT_SUBSIDENCE:
                subsidence = False
            else:
                raise ValueError(f"malformed aqueduct subsidence string {sub_str=}")

            # sea level rise percentile
            slr_perc_str = "_".join(slr_perc_list)

            # N.B. default is 95th percentile
            if slr_perc_str == "0":
                slr_percentile = 95.0
            elif slr_perc_str == "0_perc_50":
                slr_percentile = 50.0
            elif slr_perc_str == "0_perc_05":
                slr_percentile = 5.0
            else:
                raise ValueError(
                    f"malformed aqueduct sea level percentile string {slr_perc_str=}"
                )

            self.riverine = False
            self.coastal = True
            self.subsidence = subsidence
            self.model = sub_str
            self.slr_percentile = slr_percentile

        else:
            raise ValueError(
                f"do not recognise hazard {map_type=}, "
                f"{name=} must begin with either {self.RIVERINE} or {self.COASTAL}"
            )

        # attributes common to riverine and coastal maps
        self.year = int(year)
        self.return_period_years = float(return_period_years.replace("rp", ""))
        self.scenario = scenario

        return

    @property
    def without_model(self) -> str:
        """
        A name identifying the attributes of the hazard, without any model
        information (climate model / subsidence).
        """
        split_name = self.name.split("_")
        split_name.pop(2)  # index of climate model / subsidence component
        return "_".join(split_name)

    @property
    def without_RP(self) -> str:
        """
        A name identifying the attributes of the hazard, without any return
        period.

        N.B. For collapsing return periods into expected annual damages (EAD)
        it is useful to generate a name without return period information.
        """

        split_name = self.name.split("_")
        split_name.pop(4)  # index of return period element for river and coastal map names
        return "_".join(split_name)


def generate_rp_maps(names: list[str], prefix: Union[None, str] = None) -> list[ReturnPeriodMap]:
    """
    Given a list of strings, generate some ReturnPeriodMap objects. Optionally
    remove a prefix string from the input.
    """
    if prefix is not None:
        names = [re.sub(f"^{prefix}", "", name) for name in names]
    return [get_rp_map(name) for name in natural_sort(names)]


def get_rp_map(name: str) -> ReturnPeriodMap:
    """
    Given a name string, return an instance of the appropriate ReturnPeriodMap
    subclass.
    """

    # registry of implemented return period constructors
    # new return period map types must be added here
    prefix_class_map: dict[str, type[ReturnPeriodMap]] = {
        AqueductFlood.RIVERINE: AqueductFlood,
        AqueductFlood.COASTAL: AqueductFlood,
    }

    # choose constructor on name prefix
    prefix, *_ = name.split("_")

    # return a concrete subclass of ReturnPeriodMap
    return prefix_class_map[prefix](name)


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
        expected_annual_damages_path: str = snakemake.output["expected_annual_damages"]
        return_period_and_ead_path: str = snakemake.output["return_period_and_ead"]
        damage_curves_dir: str = snakemake.config["direct_damages"]["curves_dir"]
        network_type: str = snakemake.params["network_type"]
        hazard_type: str = snakemake.params["hazard_type"]
        asset_types: set[str] = set(snakemake.config["direct_damages"]["asset_types"])
    except NameError:
        raise ValueError("Must be run via snakemake.")

    OUTPUT_FILE_PATHS: tuple[str] = (damage_fraction_path, damage_cost_path, expected_annual_damages_path, return_period_and_ead_path)

    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    # load curves first so if we fail here, we've failed early
    # and we don't try and load the (potentially large) exposure file
    damage_curves = load_damage_curves(damage_curves_dir, hazard_type, asset_types)
    logging.info(f"Available damage curves: {damage_curves.keys()}")

    logging.info("Reading exposure (network/raster intersection) data")
    exposure: gpd.GeoDataFrame = gpd.read_parquet(exposure_path)
    logging.info(f"{exposure.shape=}")

    if exposure.empty:
        logging.info("No data in geometry column, writing empty files.")

        # snakemake requires that output files exist, even if empty
        for path in OUTPUT_FILE_PATHS:
            utils.write_empty_frames(path)
        sys.exit(0)  # exit gracefully so snakemake will continue

    # column groupings for data selection
    hazard_columns = [col for col in exposure.columns if col.startswith(HAZARD_PREFIX)]
    non_hazard_columns = list(set(exposure.columns) - set(hazard_columns))

    ##############################################################
    ### DAMAGE FRACTIONS (per split geometry, for all rasters) ###
    ##############################################################

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

    #######################################################################
    ### DAMAGE COST (for split, then grouped geometry, for all rasters) ###
    #######################################################################

    # multiply the damage fraction estimates by a cost to rebuild the asset
    # units are: 1 * USD/km * km = USD
    logging.info("Calculating direct damage costs")
    direct_damages_only = damage_fraction[hazard_columns] \
        .multiply(damage_fraction[REHABILITATION_COST_FIELD], axis="index") \
        .multiply(damage_fraction[SPLIT_LENGTH_FIELD], axis="index")

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

    ###################################################################################
    ### EXPECTED ANNUAL DAMAGE COST (for grouped geometry, aggregations of rasters) ###
    ###################################################################################

    # reduce climate model / subsidence to simple MIN/MAX
    # generate a mapping from a 'family' of hazards to their set of related return period maps
    model_families: dict[str, set[ReturnPeriodMap]] = defaultdict(set)
    for rp_map in generate_rp_maps(grouped_direct_damages.columns, prefix=HAZARD_PREFIX):
        model_families[rp_map.without_model].add(rp_map)  # only differ by climate model / subsidence

    logging.info("Min/max aggregating input raster models")
    for family_name, family_rp_maps in model_families.items():
        for agg_str in ("min", "max"):
            sample_map, *_ = family_rp_maps
            family_aggregation_name = sample_map.name.replace(sample_map.model, agg_str.upper())
            family_column_names: list[str] = [f"{HAZARD_PREFIX}{rp_map.name}" for rp_map in family_rp_maps]
            agg_func = getattr(grouped_direct_damages, agg_str)
            grouped_direct_damages[f"{HAZARD_PREFIX}{family_aggregation_name}"] = agg_func(axis="columns")
        grouped_direct_damages = grouped_direct_damages.drop(columns=family_column_names)

    # integrate over return periods for expected annual damages
    rp_map_families: dict[str, set[ReturnPeriodMap]] = defaultdict(set)
    for rp_map in generate_rp_maps(grouped_direct_damages.columns, prefix=HAZARD_PREFIX):
        rp_map_families[rp_map.without_RP].add(rp_map)  # only differ by return period

    expected_annual_damages = {}
    logging.info(f"Integrating {len(rp_map_families)} damage-probability curves")
    for family_name, family_rp_maps in rp_map_families.items():

        # sort by least to most probable
        sorted_rp_maps: list[ReturnPeriodMap] = sorted(family_rp_maps)

        # [0, 1] valued decimal probabilities
        probabilities: list[float] = [rp_map.annual_probability for rp_map in sorted_rp_maps]
        # family subset of grouped_direct_damages
        family_column_names: list[str] = [f"{HAZARD_PREFIX}{rp_map.name}" for rp_map in sorted_rp_maps]
        family_direct_damages: pd.DataFrame = grouped_direct_damages[family_column_names]

        # integrate the damage as a function of probability curve using Simpson's rule
        # Simpson's rule as the function to be integrated is non-linear
        expected_annual_damages[family_name] = simpson(family_direct_damages, x=probabilities, axis=1)

    #############################################
    ### JOINING, VALIDATION AND SERIALIZATION ###
    #############################################

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
    expected_annual_damages = gpd.GeoDataFrame(unsplit_subset.join(expected_annual_damages_only, validate="one_to_one", how="right"))
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

    logging.info(f"Writing out {return_period_and_ead_damages.shape=} (per split geometry, hazard RP map and hazard map (integrated RP))")
    return_period_and_ead_damages.to_parquet(return_period_and_ead_path)

    logging.info("Done calculating direct damages")
