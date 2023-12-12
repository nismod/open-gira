"""
Functions for solving the inverse geodesic problem, that is, given pairs of
points (in this case on a sphere), find the great-circle distances and (foward)
headings between them.
"""

from typing import Union

import numba
import numpy as np


@numba.njit
def forward_azimuth(
    Δλ: Union[float, np.ndarray],
    φ1: Union[float, np.ndarray],
    φ2: Union[float, np.ndarray]
):
    """
    Calculate forward azimuth as part of inverse geodesy problem (find azimuths
    and great-circle distance between pairs of points).

    Reference: https://www.movable-type.co.uk/scripts/latlong.html

    Args:
        Δλ: Difference in longitudes (radians)
        φ1: Start latitudes (radians)
        φ2: End latitudes (radians)

    Returns: Initial headings from start towards end points (radians)
    """

    y = np.sin(Δλ) * np.cos(φ2)
    x = np.cos(φ1) * np.sin(φ2) - np.sin(φ1) * np.cos(φ2) * np.cos(Δλ)

    return np.arctan2(y, x)


@numba.njit
def bearing_and_great_circle_distance(
    longitude1: Union[float, np.ndarray],
    latitude1: Union[float, np.ndarray],
    longitude2: Union[float, np.ndarray],
    latitude2: Union[float, np.ndarray]
):
    """
    Given pairs of points in decimal degrees, calculate forward azimuth in
    degrees and great-circle (spherical geometry) distance in meters.

    This implements the Haversine formula for calculating the distance, which
    is more performant but less accurate than approaches based on ellipsoids,
    such as Vincenty's algorithm. Distance errors of this function are < ~0.5%
    w.r.t. to Vincenty. Azimuth errors are similar, except in the case of
    antipodes (where the correct heading is ill-defined). This function is ~7x
    faster than pyproj.Geod.inv for ~10-1,000,000 point input arrays.

    Args:
        longitude1: Start longitudes (degrees)
        latitude1: Start latitudes (degrees)
        longitude2: End longitudes (degrees)
        latitude2: End latitudes (degrees)

    Returns:
        Forward azimuth angle (from start to end) (degrees)
        Great-circle distance between points (meters)
    """

    # convert to radians
    λ1 = np.radians(longitude1)
    φ1 = np.radians(latitude1)
    λ2 = np.radians(longitude2)
    φ2 = np.radians(latitude2)

    # latitude and longitude ranges
    Δλ = λ2 - λ1
    Δφ = φ2 - φ1

    # squares of half the chord lengths
    a = np.power(np.sin(Δφ / 2.0), 2) + np.cos(φ1) * np.cos(φ2) * np.power(np.sin(Δλ / 2.0), 2)

    # angles separating points
    θ = 2 * np.arcsin(np.sqrt(a))

    # +/- 10km depending on (mostly) latitude
    earth_radius_m = 6371000

    # great circle distances, l = rθ
    distance_m = earth_radius_m * θ

    return np.rad2deg(forward_azimuth(Δλ, φ1, φ2)), distance_m
