"""
Given an exposure estimate and some damage curves, calculate the damage
fraction for exposed assets.
"""

from abc import ABC, abstractmethod
from collections import defaultdict
import logging
from typing import Union
import warnings
import re

import geopandas as gpd
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.interpolate import interp1d
from scipy.integrate import simpson
from tqdm import tqdm

def natural_sort(to_sort):
    return sorted(to_sort, key=lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)])


# exposure table hazard intensity fields expected to be prefixed as such
HAZARD_PREFIX = "hazard"
# exposure table field containing the cost to rebuild per unit length
REHABILITATION_COST = 200_000
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
    def without_rp(self):
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


class StormWindspeed(ReturnPeriodMap):
    """
    Class holding information about STORM return period windspeed maps.

    Each point in these raster maps is a wind speed for a given
    combination of e.g. scenaro, climate model, year, return period
    (probability).
    """

    STORM_GCMS = ("CMCC-CM2-VHR4", "CNRM-CM6-1-HR", "EC-Earth3P-HR", "HadGEM3-GC31-HM")

    def __init__(self, name: str):
        """
        Infer attributes from name.

        Arguments:
            name (str): Name string expected to be in one of the following formats:
                Current:
                    hazard-storm_{year}_{rcp}_rp{storm_rp}_{min/max/model}
        """
        slug = name.replace("hazard-storm_", "")

        # store the original string for later use
        self.name = name

        try:
            year, rcp, rp, gcm = slug.split("_")
        except ValueError as ex:
            logging.info(f"{slug=}")
            raise ex

        self.return_period_years = int(rp.replace("rp", ""))
        self.model = gcm
        self.year = year
        self.scenario = rcp

    @classmethod
    def from_raw(cls, name: str):
        """
        Infer attributes from name.

        Arguments:
            name (str): Name string expected to be in one of the following formats:
                Current:
                    STORM_FIXED_RETURN_PERIODS_constant_{storm_rp}_YR_RP
                Future:
                    STORM_FIXED_RETURN_PERIODS_{storm_gcm}_{storm_rp}_YR_RP
        """
        slug = name.replace("STORM_FIXED_RETURN_PERIODS_", "").replace("_YR_RP", "")
        gcm, rp = slug.split("_")
        if gcm == "constant":
            year = 2020
            rcp = "baseline"
        else:
            year = 2050
            rcp = "rcp8p5"
        return StormWindspeed(f"hazard-storm_{year}_{rcp}_rp{rp}_{gcm}")

    @property
    def without_model(self) -> str:
        """
        A name identifying the attributes of the hazard, without any model
        information (climate model / subsidence).
        """
        return f"hazard-storm_{self.year}_{self.scenario}_rp{self.return_period_years}"

    @property
    def without_rp(self) -> str:
        """
        A name identifying the attributes of the hazard, without any return
        period.

        N.B. For collapsing return periods into expected annual damages (EAD)
        it is useful to generate a name without return period information.
        """
        return f"hazard-storm_{self.year}_{self.scenario}_{self.model}_ead"


def get_rp_map(name: str) -> ReturnPeriodMap:
    """
    Given a name string, return an instance of the appropriate ReturnPeriodMap
    subclass.
    """
    # return a concrete subclass of ReturnPeriodMap
    return StormWindspeed(name)

def generate_rp_maps(names: list[str], prefix: Union[None, str] = None) -> list[ReturnPeriodMap]:
    """
    Given a list of strings, generate some ReturnPeriodMap objects. Optionally
    remove a prefix string from the input.
    """
    if prefix is not None:
        names = [re.sub(f"^{prefix}", "", name) for name in names]
    return [get_rp_map(name) for name in natural_sort(names)]

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
    windspeed_damage = pd.DataFrame({
        "windspeed_kmh":   [ 0.0, 120, 150, 180, 200, 220, 230, 250],
        "damage_fraction": [ 0.0, 0.0, 0.1, 0.2, 0.4, 0.8, 0.9, 1.0],
    })
    windspeed_damage["windspeed_ms"] = windspeed_damage.windspeed_kmh * (1000 / 3600)
    windspeed_damage = windspeed_damage[["windspeed_ms", "damage_fraction"]]
    return {
        "transmission": windspeed_damage
    }


def calculate_damages(exposure: pd.DataFrame) -> tuple[pd.DataFrame]:
    # column groupings for data selection
    hazard_columns = [col for col in exposure.columns if col.startswith(HAZARD_PREFIX)]
    non_hazard_columns = list(set(exposure.columns) - set(hazard_columns))
    asset_exposure = exposure.loc[:, hazard_columns]

    # apply damage_curve function to exposure table
    # the return value of interpolated_damage_curve is a numpy array
    damage_fraction_data = pd.DataFrame(
        interpolated_damage_curve(asset_exposure),
        index=asset_exposure.index,
        columns=asset_exposure.columns
    )
    # store the computed direct damages and any columns we started with
    damage_fraction = pd.concat(
        [
            damage_fraction_data,
            exposure.loc[:, non_hazard_columns]
        ],
        axis="columns"
    )

    # multiply the damage fraction estimates by a cost to rebuild the asset
    # units are: 1 * USD/km * km = USD
    direct_damages_data = damage_fraction_data \
        .multiply(REHABILITATION_COST, axis="index") \
        .multiply(damage_fraction[SPLIT_LENGTH_FIELD], axis="index")

    # join the other fields with the direct damage estimates
    direct_damages = pd.concat(
        [direct_damages_data, damage_fraction["edge_id"]],
        axis="columns"
    ).set_index("edge_id")

    # grouping on edge_id, sum all direct damage estimates to give a total dollar cost per edge
    direct_damages = direct_damages.groupby(direct_damages.index).sum()


    # reduce climate model to simple MIN/MEAN/MAX
    # generate a mapping from a 'family' of hazards to their set of related return period maps
    model_families: dict[str, set[ReturnPeriodMap]] = defaultdict(set)
    for rp_map in generate_rp_maps(direct_damages.columns):
        model_families[rp_map.without_model].add(rp_map)  # only differ by climate model / subsidence

    for family_name, family_rp_maps in model_families.items():
        for agg_str in ("min", "mean", "max"):
            family_aggregation_name = f"{family_name}_{agg_str}"
            family_column_names: list[str] = [rp_map.name for rp_map in family_rp_maps]
            if agg_str == "min":
                agg = direct_damages[family_column_names].min(axis="columns")
            elif agg_str == "max":
                agg = direct_damages[family_column_names].max(axis="columns")
            elif agg_str == "mean":
                agg = direct_damages[family_column_names].mean(axis="columns")
            else:
                assert False, "unrecognised aggregation"
            direct_damages[family_aggregation_name] = agg
        direct_damages = direct_damages.drop(columns=family_column_names)


    # integrate over return periods for expected annual damages
    rp_map_families: dict[str, set[ReturnPeriodMap]] = defaultdict(set)
    for rp_map in generate_rp_maps(direct_damages.columns):
        rp_map_families[rp_map.without_rp].add(rp_map)  # only differ by return period

    expected_annual_damages = {}
    for family_name, rp_maps in rp_map_families.items():
        sorted_rp_maps: list[ReturnPeriodMap] = sorted(rp_maps)

        # [0, 1] valued decimal probabilities (least to most probable now we've sorted)
        probabilities: list[float] = [rp_map.annual_probability for rp_map in sorted_rp_maps]
        # family subset of direct_damages
        family_column_names: list[str] = [rp_map.name for rp_map in sorted_rp_maps]
        family_direct_damages: pd.DataFrame = direct_damages[family_column_names]

        # integrate the damage as a function of probability curve using Simpson's rule
        # Simpson's rule as the function to be integrated is non-linear
        expected_annual_damages[family_name] = simpson(family_direct_damages, x=probabilities, axis=1)

    expected_annual_damages = pd.DataFrame(data=expected_annual_damages, index=direct_damages.index)
    return direct_damages, expected_annual_damages

if __name__ == "__main__":

    try:
        unsplit_path: str = "results/grid.geoparquet"
        exposure_path: str = "results/grid_split.geoparquet"
        damage_path: str = "results/grid_damage.geoparquet"
        damage_curves_dir: str = ""
        network_type: str = "transmission"
        hazard_type: str = "storm"
        asset_types: set[str] = "transmission"
    except NameError:
        raise ValueError("Must be run via snakemake.")

    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")

    # load curves first so if we fail here, we've failed early
    # and we don't try and load the (potentially large) exposure file
    damage_curves = load_damage_curves(damage_curves_dir, hazard_type, asset_types)

    # calculate damages
    damage_curve: pd.DataFrame = damage_curves["transmission"]
    print(damage_curve)

    # create interpolated damage curve for given asset type
    hazard_intensity, damage_fraction = damage_curve.iloc[:, 0], damage_curve.iloc[:, 1]
    # if curve of length n, where x < x_0, y = y_0 and where x > x_n, y = y_n
    bounds = tuple(f(damage_fraction) for f in (min, max))
    interpolated_damage_curve = interp1d(hazard_intensity, damage_fraction, kind='linear', fill_value=bounds, bounds_error=False)


    exp_pf = pq.ParquetFile(exposure_path)
    batch_size = 20_000

    for i, exp_batch in tqdm(enumerate(exp_pf.iter_batches(batch_size))):
        exposure = exp_batch.to_pandas()
        cols = []
        for c in exposure.columns:
            if "STORM" in c:
                c = StormWindspeed.from_raw(c.replace(".tif", "")).name
            cols.append(c)
        exposure.columns = cols

        dd, ead = calculate_damages(exposure)
        dd.join(ead).to_parquet(damage_path.replace(".geo", f"_{i}."))



    logging.info("Done")
