"""
Functions for reconstructing wind fields from track data.
"""

import numpy as np
import pyproj

import geopandas as gpd
import pandas as pd


def power_law_scale_factors(z0: np.ndarray, z1: float, z2: float) -> np.ndarray:
    """
    Calculate factors to scale wind speeds from one height to another, given
    surface roughness information.

    Wieringa 1992, Journal of Wind Engineering and Industrial Aerodynamics
    Volume 41, Issues 1–3, October 1992, Pages 357-368, Equation 3

    This approach is rudimentary, but requires few inputs. Best for scaling
    from top of boundary layer (~2km) to a few hundred metres above ground.
    Logarithmic scaling better suited for 100m -> ground level.

    Args:
        z0: Surface roughness lengths in metres
        z1: Height of desired wind speeds above ground in metres
        z2: Height of wind speeds, v2, above ground in metres

    Returns:
        Scale factors, f, for v2 = f * v1
    """

    p: np.ndarray = 1 / np.log(np.sqrt(z1 * z2) / z0)

    return np.power(z1 / z2, p)


def holland_wind_model(
    RMW_m: float,
    V_max_ms: float,
    p_pa: float,
    p_env_pa: float,
    r_m: np.ndarray,
    phi_deg: float,
) -> np.ndarray:
    """
    Calculate gradient-level wind speed at points some distance from a cyclone eye location.

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

    # case where Delta_P is zero will raise ZeroDivisionError
    Delta_P = p_env_pa - p_pa

    # beta parameter, how sharp the V(r) profile around the eye wall is
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


def advective_vector(
    eye_heading_deg: float,
    eye_speed_ms: float,
    hemisphere: int,
    alpha: float = 0.56,
    beta: float = 19.2,
) -> np.complex128:
    """
    Calculate the advective wind vector

    For reconstruction of realistic advective wind component, (and rationale
    for alpha and beta) see section 2 of: Lin, N., and D. Chavas (2012), On
    hurricane parametric wind and applications in storm surge modeling, J.
    Geophys. Res., 117, D09120, doi:10.1029/2011JD017126

    Arguments:
        eye_heading_deg (float): Heading of eye in degrees clockwise from north
        eye_speed_ms (float): Speed of eye in metres per second
        hemisphere (int): +1 for northern, -1 for southern
        alpha (float): Fractional reduction of advective wind speed from eye speed
        beta (float): Degrees advective winds tend to rotate past storm track (in direction of rotation)

    Returns:
        np.complex128: Advective wind vector
    """

    # bearing of advective component (storm track heading with beta correction)
    phi_a: float = np.radians(eye_heading_deg - hemisphere * beta)

    # absolute magnitude of vector is eye speed decreased by alpha factor
    mag_v_a: float = eye_speed_ms * alpha

    # find components
    return mag_v_a * np.sin(phi_a) + mag_v_a * np.cos(phi_a) * 1j


def rotational_field(
    longitude: np.ndarray,
    latitude: np.ndarray,
    eye_long: float,
    eye_lat: float,
    radius_to_max_winds_m: float,
    max_wind_speed_ms: float,
    min_pressure_pa: float,
    env_pressure_pa: float,
) -> np.ndarray:
    """
    Calculate the rotational component of a storm's vector wind field

    Args:
        longitude (np.ndarray[float]): Grid values to evaluate on
        latitude (np.ndarray[float]): Grid values to evaluate on
        eye_long (float): Location of eye in degrees
        eye_lat (float): Location of eye in degrees
        radius_to_max_winds_m (float): Distance from eye centre to maximum wind speed in metres
        max_wind_speed_ms (float): Maximum linear wind speed (relative to storm eye)
        min_pressure_pa (float): Minimum pressure in storm eye in Pascals
        env_pressure_pa (float): Environmental pressure, typical for this locale, in Pascals

    Returns:
        np.ndarray[complex]: Grid of wind vectors
    """

    X, Y = np.meshgrid(longitude, latitude)
    grid_shape = X.shape  # or Y.shape

    # forward azimuth angle and distances from grid points to track eye
    geod_wgs84: pyproj.Geod = pyproj.CRS("epsg:4326").get_geod()
    grid_to_eye_azimuth_deg, _, radius_m = geod_wgs84.inv(
        X.ravel(),
        Y.ravel(),
        np.full(len(X.ravel()), eye_long),
        np.full(len(Y.ravel()), eye_lat),
    )

    # magnitude of rotational wind component
    mag_v_r: np.ndarray = holland_wind_model(
        radius_to_max_winds_m,
        max_wind_speed_ms,
        min_pressure_pa,
        env_pressure_pa,
        radius_m.reshape(grid_shape),
        eye_lat
    )

    # azimuth of rotational component is tangent to radius, with direction set by hemisphere
    phi_r: np.ndarray = np.radians(grid_to_eye_azimuth_deg.reshape(grid_shape) + np.sign(eye_lat) * 90)

    # find components of vector at each pixel
    return mag_v_r * np.sin(phi_r) + mag_v_r * np.cos(phi_r) * 1j


def interpolate_track(track: gpd.GeoDataFrame, frequency: str = "1H") -> gpd.GeoDataFrame:
    """
    Interpolate storm track data.

    Arguments:
        track (gpd.GeoDataFrame): Storm track with at least the following
            columns: geometry, min_pressure_hpa, max_wind_speed_ms,
            radius_to_max_winds_km, timestep. Must have a DatetimeIndex.
        frequency (str): If given track with DatetimeIndex, interpolate to
            resolution given by this pandas frequency string

    Returns:
        gpd.GeoDataFrame: Track with min_pressure_hpa, max_wind_speed_ms,
            radius_to_max_winds_km and geometry columns interpolated.
    """

    if len(track) == 0:
        raise ValueError("No track data")
    elif len(track) == 1:
        # not enough data to interpolate between, short circuit
        return track
    elif len(track) == 2:
        interp_method = "linear"
    else:
        interp_method = "quadratic"

    if isinstance(track.index, pd.DatetimeIndex):
        interp_index = pd.date_range(track.index[0], track.index[-1], freq=frequency)
    else:
        raise ValueError("tracks must have a datetime index to interpolate")

    track["x"] = track.geometry.x
    track["y"] = track.geometry.y
    track = track.drop(columns="geometry").copy()

    interp_cols = [
        "min_pressure_hpa",
        "max_wind_speed_ms",
        "radius_to_max_winds_km",
        "x",
        "y",
    ]

    interp_domain = pd.DataFrame(
        index=interp_index,
        data=np.full((len(interp_index), len(track.columns)), np.nan),
        columns=track.columns,
    )

    # merge dataframes, on index collision, keep the non-NaN values from track
    interp_track = track.combine_first(interp_domain).sort_index()

    # interpolate over numeric value of index
    interp_track.loc[:, interp_cols] = interp_track.loc[:, interp_cols].interpolate(method=interp_method)

    interp_track["geometry"] = gpd.points_from_xy(
        interp_track.x, interp_track.y, crs="EPSG:4326"
    )

    return gpd.GeoDataFrame(interp_track).drop(columns=["x", "y"])
