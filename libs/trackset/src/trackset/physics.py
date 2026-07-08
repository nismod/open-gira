"""
Parametric models for tropical cyclone structure.

Some synthetic track sets (e.g. CHAZ) report only position and maximum wind
speed; the radius to maximum winds and minimum pressure required by wind
field models must be inferred. These fits are ported from
https://github.com/thomas-fred/chaz-track-parser.
"""

import numpy as np

# Environmental pressure values in hPa / mbar (standard estimate of background
# pressure away from the cyclone) are taken from the AIR hurricane model,
# table 3 in Butke (2012). Available at:
# https://www.air-worldwide.com/publications/air-currents/2012/
# the-pressures-on-increased-realism-in-tropical-cyclone-wind-speeds-through-attention-to-environmental-pressure/
ENV_PRESSURE = {
    "NI": 1006.5,
    "SA": 1014.1,
    "NA": 1014.1,
    "EP": 1008.8,
    "SI": 1010.6,
    "SP": 1008.1,
    "WP": 1008.3,
}

# Rotation speed of the Earth in rad/s
OMEGA = (2 * np.pi) / (24 * 60 * 60)

HECTOPASCALS_PER_PASCAL = 0.01

# See Emanuel 1986 for a discussion of minimum possible pressures (circa
# 850 hPa); climate change may reduce this slightly, but very unlikely to
# pass 800 hPa. https://doi.org/10.1175/1520-0469(1986)043%3C0585:AASITF%3E2.0.CO;2
MIN_PLAUSIBLE_PRESSURE_HPA = 800.0
MAX_PLAUSIBLE_PRESSURE_HPA = 1050.0


def r_max_willoughby_2004(v_max: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """
    Tropical cyclone radius to maximum wind (km) as a function of maximum
    sustained wind speed (m/s) and latitude (degrees).

    Log-linear fit to aircraft observation data.

    Source: https://doi.org/10.1175/MWR2831.1, equation 12.1.

    N.B. The source does not include taking the absolute value of the
    latitude, but it cannot be correct for both hemispheres otherwise.
    """

    return 46.29 * np.exp(-0.0153 * v_max + 0.0166 * np.abs(phi))


def coriolis(phi: np.ndarray) -> np.ndarray:
    """
    Coriolis parameter (radians / second) for given latitudes (degrees).
    """

    return np.abs(2 * OMEGA * np.sin(np.radians(phi)))


def b_vickery_wadhera_2008(phi: np.ndarray, r_max: np.ndarray) -> np.ndarray:
    """
    Holland's B (shape) parameter as a function of latitude (degrees) and
    radius to maximum winds (metres), courtesy of a fit in Vickery and
    Wadhera 2008.

    Source: https://doi.org/10.1175/2008JAMC1837.1, equation 23.
    """

    return 1.833 - 0.326 * np.sqrt(coriolis(phi) * r_max)


def p_min_holland_1980(
    p_env: float | np.ndarray,
    v_max: np.ndarray,
    r_max: np.ndarray,
    phi: np.ndarray,
    rho: float | np.ndarray = 1.15,
) -> np.ndarray:
    """
    Infer minimum pressure of a tropical cyclone from basic information.

    When possible, it is best to use a model with a secondary shape
    parameter, such as radius to 34kt winds (see Chavas 2024). In the absence
    of this additional shape data, this function offers a rudimentary
    estimate of pressure.

    The derivation of the implemented equation is as follows:

    Gradient wind balance (see http://noaa-ocs-modeling.github.io/PaHM/html/models.html):
    v(r)**2 + frv(r) = r/rho * dp(r)/dr, (1)

    Holland 1980 pressure profile:
    p(r) = p_min + (p_env - p_min) * e**(-(R_max/r)**B), (2)

    Find the derivative of (2) where r = r_max:
    dp(r)/dr = B(p_env - p_min) / (e * r_max), (3)

    Then, substitute (3) into (1), again with r = r_max and rearrange to give:
    p_min = p_env - (rho * e * (v_max**2 + f * r_max * v_max)) / B

    Args:
        p_env: Environmental pressure (well away from influence of TC), hectopascals
        v_max: Maximum rotational wind speed, metres / second
        r_max: Radial distance from eye to maximum wind speed, metres
        phi: Latitude of eye locations, degrees
        rho: Density of air, kilograms per cubic metre

    Returns:
        An estimate of minimum pressure, hectopascals. Implausible values
        (outside 800 - 1050 hPa) are returned as NaN.
    """

    lhs_wind_balance = v_max**2 + coriolis(phi) * r_max * v_max
    numerator = rho * np.e * lhs_wind_balance
    b = b_vickery_wadhera_2008(phi, r_max)
    p_min = p_env - HECTOPASCALS_PER_PASCAL * (numerator / b)
    return np.where(
        (p_min < MIN_PLAUSIBLE_PRESSURE_HPA) | (p_min > MAX_PLAUSIBLE_PRESSURE_HPA),
        np.nan,
        p_min,
    )
