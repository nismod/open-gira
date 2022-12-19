"""
Functions for reconstructing wind fields from track data.
"""

import numpy as np


def holland_wind_model(
    RMW_m: float,
    V_max_ms: float,
    p_pa: float,
    p_env_pa: float,
    r_m: np.ndarray,
    phi_deg: float,
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
        V_max_ms (float): Maximum surface wind speed (minus background flow) in
            meters per second
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

    Omega = (2 * np.pi) / (24 * 60 * 60)  # rotation speed of the Earth in rad/s
    f = np.abs(2 * Omega * np.sin(np.radians(phi_deg)))  # Coriolis parameter

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
