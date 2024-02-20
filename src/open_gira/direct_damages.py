"""
Functionality to assist calculating the direct damages to assets due to hazards.
"""

from abc import ABC, abstractmethod
import logging
import re
from typing import Union

import geopandas as gpd
import pandas as pd
from scipy.interpolate import interp1d

from open_gira import fields
from open_gira.utils import natural_sort


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


def rail_rehab_cost(row: pd.Series, rehab_cost: pd.DataFrame) -> float:
    """
    Determine the cost of rehabilitation for a given rail segment (row).

    Args:
        row: Railway segment
        rehab_cost: Table of rehabilitation costs for various rail asset types

    Returns:
        Cost in units of USD per km.
    """

    data = rehab_cost.loc[rehab_cost.asset_type == row.asset_type, "rehab_cost_USD_per_km"]
    if data.empty:
        raise ValueError(f"Missing rehabilitation cost data for {row.asset_type}, please amend.")

    # unpack single element series
    cost, = data

    return float(cost)


def road_rehab_cost(row: pd.Series, rehab_cost: pd.DataFrame) -> float:
    """
    Fetch the cost of rehabilitation for a given road segment (row).

    Args:
        row: Road segment
        rehab_cost: Table of rehabilitation costs for various road asset types

    Returns:
        Cost in units of USD per km per lane.
    """

    data = rehab_cost.loc[rehab_cost.asset_type == row.asset_type, "rehab_cost_USD_per_km_per_lane"]
    if data.empty:
        raise ValueError(f"Missing rehabilitation cost data for {row.asset_type}, please amend.")

    # unpack single element series
    cost, = data

    return float(cost)


def annotate_rehab_cost(edges: gpd.GeoDataFrame, network_type: str, rehab_cost: pd.DataFrame) -> gpd.GeoDataFrame:
    """
    Label edges with rehabilitation costs to be used in subsequent direct
    damage cost estimate.

    Args:
        edges: Edges to label.
        network_type: Category string, currently either road or rail.
        rehab_cost: Table containing rehabilitation cost data per km. See lookup
            functions for more requirements.

    Returns:
        Edges labelled with rehabilitation cost data.
    """

    if network_type == "road":
        edges[fields.REHAB_COST] = edges.apply(road_rehab_cost, axis=1, args=(rehab_cost,)) * edges.lanes
    elif network_type == "rail":
        edges[fields.REHAB_COST] = edges.apply(rail_rehab_cost, axis=1, args=(rehab_cost,))
    else:
        raise ValueError(f"No lookup function available for {network_type=}")

    return edges


def direct_damage(
    exposure: gpd.GeoDataFrame,
    damage_curves: dict[str, pd.DataFrame],
    hazard_columns: list[str],
    non_hazard_columns: list[str],
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Calculate direct damages for exposed edge assets. Take hazard intensity,
    lookup damage fraction from a damage curve, multiply by a rehabilitation
    cost.

    Args:
        exposure: Table containing exposed assets, likely edges split on a
            raster grid (i.e. `edge_id` is not necessarily unique).
        damage_curves: Relationship between hazard intensity and damage
            fraction, keyed by `asset_type`.
        hazard_columns: Columns in `exposure` which denote hazard intensities.
        non_hazard_columns: Columns in `exposure` which do not denote hazard
            intensities.

    Returns:
        Direct damage fraction, rows are splits of edges.
        Direct damage cost, rows are edges and `edge_id` should now be unique.
    """

    ##########################################################
    # DAMAGE FRACTIONS (per split geometry, for all rasters) #
    ##########################################################

    # calculate damages for assets we have damage curves for
    damage_fraction_by_asset_type = []
    logging.info(f"Exposed assets {natural_sort(set(exposure.asset_type))}")
    for asset_type in natural_sort(set(exposure.asset_type) & set(damage_curves.keys())):

        damage_curve: pd.DataFrame = damage_curves[asset_type]

        # pick out rows of asset type and columns of hazard intensity
        asset_type_mask: gpd.GeoDataFrame = exposure.asset_type == asset_type
        asset_exposure: pd.DataFrame = pd.DataFrame(exposure.loc[asset_type_mask, hazard_columns])

        logging.info(f"Processing {asset_type=}, with n={len(asset_exposure)}")

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

    # join the other fields with the direct damage estimates
    logging.info("Unifying rasterised segments and summing damage costs")

    # grouping on edge_id, sum all direct damage estimates to give a total dollar cost per edge
    direct_damages = pd.concat(
        [direct_damages_only, damage_fraction["edge_id"]],
        axis="columns"
    ).set_index("edge_id")
    grouped_direct_damages = direct_damages.groupby(direct_damages.index).sum()

    return damage_fraction, grouped_direct_damages