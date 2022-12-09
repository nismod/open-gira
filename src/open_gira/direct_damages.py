"""
Functionality to assist calculating the direct damages to assets due to hazards.
"""

from abc import ABC, abstractmethod
import re
from typing import Union

import numpy as np

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


def holland_wind_model(
    RMW_m: float,
    V_max_ms: float,
    p_pa: float,
    p_env_pa: float,
    r_m: np.ndarray,
    psi_deg: float,
) -> np.ndarray:
    """
    Calculate wind speed at points some distance from a cyclone eye location.

    References
    ----------
    - Lin and Chavas (2012)
      https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2011JD017126
    - Holland (1980)
      https://doi.org/10.1175/1520-0493(1980)108%3C1212:AAMOTW%3E2.0.CO;2

    See in particular Section 3.2 in Lin and Chavas (2012).

    Args:
        RMW_m (float): Radius to max wind speeds in meters
        V_max_ms (float): Maximum wind speed in meters per second
        p_pa (float): Pressure of eye in Pascals
        p_env_pa (float): 'Background' atmospheric pressure in Pascals
        r_m (np.ndarray): Radii in meters to calculate wind speeds for
        phi_deg (float): Latitude of eye in degrees

    Returns:
        np.ndarray: Wind speeds given other params for each radii in input r_m
    """

    M = 0.02897  # molar mass of (dry) air, kg/mol
    R = 8.314  # gas constant, J/K/mol
    T = 293  # temperature estimate, K
    rho = (p_pa * M) / (R * T)  # kg/m^3

    Omega = 7.292E-5
    f = np.abs(2 * Omega * np.sin(np.radians(psi_deg)))  # Coriolis parameter

    # case where (pressure_env_hpa == pressure_hpa) so Delta_P is zero will raise ZeroDivisionError
    Delta_P = p_env_pa - p_pa

    B = (
        np.power(V_max_ms, 2) * np.e * rho
        + f * V_max_ms * RMW_m * np.e * rho
    ) / Delta_P

    V = (
        np.sqrt(
            (
                # case where r_m is zero will raise ZeroDivisionError
                (
                    np.power(RMW_m / r_m, B)
                    * B
                    * Delta_P
                    * np.exp(0 - (RMW_m / r_m) ** B)
                )
                + (np.power(r_m, 2) * np.power(f, 2) / 4)
            )
            / rho
        )
        - (f * r_m) / 2
    )

    return np.clip(V, 0, None)  # clip negative values to zero


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
