"""
For reconstructing wind fields from track data.
"""

import numba
import numpy as np
import xarray as xr

import geopandas as gpd
import pandas as pd

from open_gira.geodesic import bearing_and_great_circle_distance


# dimensions on which we have maximum wind speeds for
WIND_COORDS: dict[str, type] = {
    "event_id": str,
    "latitude": float,
    "longitude": float,
}

# Environmental pressure values in hPa / mbar (standard estimate of background
# pressure away from the cyclone) are taken from the AIR hurricane model, table
# 3 in Butke (2012). Available at:
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


def empty_wind_da() -> xr.DataArray:
    """
    Return a maxium wind field dataarray with a schema but no data.

    N.B. To concatenate xarray objects, they must all share the same
    coordinate variables.
    """
    da = xr.DataArray(
        data=np.full((0,) * len(WIND_COORDS), np.nan),
        coords={
            name: np.array([], dtype=dtype)
            for name, dtype in WIND_COORDS.items()
        },
        name="max_wind_speed"
    )
    return da


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


@numba.njit(error_model="numpy")  # allow divide by zero
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


@numba.njit
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
        eye_heading_deg: Heading of eye in degrees clockwise from north
        eye_speed_ms: Speed of eye in metres per second
        hemisphere: +1 for northern, -1 for southern
        alpha: Fractional reduction of advective wind speed from eye speed
        beta: Degrees advective winds tend to rotate past storm track (in direction of rotation)

    Returns:
        np.complex128: Advective wind vector
    """

    # bearing of advective component (storm track heading with beta correction)
    phi_a: float = np.radians(eye_heading_deg - hemisphere * beta)

    # absolute magnitude of vector is eye speed decreased by alpha factor
    mag_v_a: float = eye_speed_ms * alpha

    # find components
    return mag_v_a * np.sin(phi_a) + mag_v_a * np.cos(phi_a) * 1j


@numba.njit
def sigmoid_decay(x: np.ndarray, midpoint: float, slope: float) -> np.ndarray:
    """
    Transform input by a sigmoid shape decay to return values in range [0, 1].

    Args:
        x: Input array to transform
        midpoint: Decay midpoint in x
        slope: Larger values decay faster
    """
    return 0.5 * (1 + np.tanh(slope * (midpoint - x)))


def estimate_wind_field(
    longitude: np.ndarray,
    latitude: np.ndarray,
    eye_long: float,
    eye_lat: float,
    radius_to_max_winds_m: float,
    max_wind_speed_ms: float,
    min_pressure_pa: float,
    env_pressure_pa: float,
    advection_azimuth_deg: float,
    eye_speed_ms: float,
) -> np.ndarray:
    """
    Given a spatial domain and tropical cyclone attributes, estimate a vector wind field.

    The rotational component uses a modified Holland model. The advective
    component (from the storm's translational motion) is modelled to fall off
    from approximately 10 maximum wind radii. These two components are vector
    added and returned.

    Args:
        longitude: Grid values to evaluate on
        latitude: Grid values to evaluate on
        eye_long: Location of eye in degrees
        eye_lat: Location of eye in degrees
        radius_to_max_winds_m: Distance from eye centre to maximum wind speed in metres
        max_wind_speed_ms: Maximum wind speed (relative to ground)
        min_pressure_pa: Minimum pressure in storm eye in Pascals
        env_pressure_pa: Environmental pressure, typical for this locale, in Pascals
        eye_heading_deg: Heading of eye in degrees clockwise from north
        eye_speed_ms: Speed of eye in metres per second

    Returns:
        Grid of wind vectors
    """

    # check inputs
    assert 0 < max_wind_speed_ms < 130
    assert 0 < radius_to_max_winds_m < 1500000
    assert 75000 < min_pressure_pa < 102000
    assert 0 <= eye_speed_ms < 40

    # clip eye speed to a maximum of 30ms-1
    # greater than this is non-physical, and probably the result of a data error
    # we do not want to propagate such an error to our advective wind field
    adv_vector: np.complex128 = advective_vector(
        advection_azimuth_deg,
        eye_speed_ms,
        np.sign(eye_lat),
    )

    # maximum wind speed, less advective component
    # this is the maximum tangential wind speed in the eye's non-rotating reference frame
    max_wind_speed_relative_to_eye_ms: float = max_wind_speed_ms - np.abs(adv_vector)

    X, Y = np.meshgrid(longitude, latitude)
    grid_shape = X.shape  # or Y.shape

    # forward azimuth angle and distances from grid points to track eye
    grid_to_eye_azimuth_deg, radius_m = bearing_and_great_circle_distance(
        X.ravel(),
        Y.ravel(),
        np.full(len(X.ravel()), eye_long),
        np.full(len(Y.ravel()), eye_lat),
    )

    distance_to_eye_grid_m = radius_m.reshape(grid_shape)

    # decay effect of advective field from maximum at storm eye out to zero at 1,000km radius
    # we shouldn't claim any authority on winds outside the vicinity of the storm
    adv_field: np.ndarray = adv_vector * sigmoid_decay(distance_to_eye_grid_m / 1_000, 500, 0.004)

    # magnitude of rotational wind component
    mag_v_r: np.ndarray = holland_wind_model(
        radius_to_max_winds_m,
        max_wind_speed_relative_to_eye_ms,
        min_pressure_pa,
        env_pressure_pa,
        distance_to_eye_grid_m,
        eye_lat
    )

    # azimuth of rotational component is tangent to radius, with direction set by hemisphere
    phi_r: np.ndarray = np.radians(grid_to_eye_azimuth_deg.reshape(grid_shape) + np.sign(eye_lat) * 90)

    # find components of vector at each pixel
    rot_field = mag_v_r * np.sin(phi_r) + mag_v_r * np.cos(phi_r) * 1j

    return adv_field + rot_field


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

    # don't interpolate over some sensible duration
    # note that we limit by an integer number of timesteps, not a time
    # so our implicit assumption is that the input index is equally spaced
    max_steps_to_fill: int = np.round(pd.Timedelta("6H") / pd.Timedelta(frequency)).astype(int)

    # interpolate over numeric value of index
    interp_track.loc[:, interp_cols] = interp_track.loc[:, interp_cols].interpolate(
        method=interp_method,
        limit=max_steps_to_fill
    )

    # fail if we still have NaN values (probably the time gaps exceed `max_steps_to_fill`
    assert not interp_track.loc[:, interp_cols].isnull().values.any()

    interp_track["geometry"] = gpd.points_from_xy(
        interp_track.x, interp_track.y, crs="EPSG:4326"
    )

    return gpd.GeoDataFrame(interp_track).drop(columns=["x", "y"])
