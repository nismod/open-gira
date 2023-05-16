"""
Process storm data and write out maximum wind speeds at each grid pixel for
each storm.
"""

import os
import multiprocessing
import logging
from typing import Optional
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
import rioxarray
import xarray as xr

from open_gira.wind import (
    advective_vector, rotational_field, interpolate_track,
    power_law_scale_factors, empty_wind_da, WIND_COORDS, ENV_PRESSURE
)
from plot_wind_fields import plot_contours, animate_track, plot_downscale_factors


logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

# wind speed altitudes
GRADIENT_LEVEL_METRES = 20
SURFACE_LEVEL_METRES = 10


def cleanup(output_path: str, downscale_factors_plot_path: str):
    """
    If we don't have a network, or tracks and we can't continue, write empty
    files and quit.
    """
    empty_wind_da().to_netcdf(output_path)
    os.system(f"touch {downscale_factors_plot_path}")
    sys.exit(0)


def process_track(
    track: pd.core.groupby.generic.DataFrameGroupBy,
    longitude: np.ndarray,
    latitude: np.ndarray,
    downscaling_factors: np.ndarray,
    plot_max_wind: bool,
    plot_animation: bool,
    plot_dir: Optional[str]
) -> tuple[str, np.ndarray]:
    """
    Interpolate a track, reconstruct the advective and rotational vector wind
    fields, sum them and take the maximum of the wind vector magnitude across
    time. Optionally plot the wind fields and save to disk.

    Args:
        track: Subset of DataFrame describing a track. Must have a temporal
            index and the following fields: `min_pressure_hpa`,
            `max_wind_speed_ms`, `radius_to_max_winds_km`.
        longitude: Longitude values to construct evaluation grid
        latitude: Latitude values to construct evaluation grid
        downscaling_factors: Factors to bring gradient-level winds to surface.
        plot_max_wind: Whether to plot max wind fields
        plot_animation: Whether to plot wind field evolution
        plot_dir: Where to save optional plots.

    Returns:
        str: Track identifier
        np.ndarray: 2D array of maximum wind speed experienced at each grid pixel
    """

    track_id, = set(track.track_id)

    logging.info(track_id)

    grid_shape: tuple[int, int] = (len(latitude), len(longitude))

    # we can't calculate the advective component without at least two points
    if len(track) == 1:
        return track_id, np.zeros(grid_shape)

    # basin of first record for storm track (storm genesis for synthetic tracks)
    basin: str = track.iloc[0, track.columns.get_loc("basin_id")]

    # interpolate track (avoid 'doughnut effect' of wind field from infrequent eye observations)
    track: gpd.GeoDataFrame = interpolate_track(track)

    geod_wgs84: pyproj.Geod = pyproj.CRS("epsg:4326").get_geod()

    # forward azimuth angle and distances from track eye to next track eye
    advection_azimuth_deg, _, eye_step_distance_m = geod_wgs84.inv(
        track.geometry.x.iloc[:-1],
        track.geometry.y.iloc[:-1],
        track.geometry.x.iloc[1:],
        track.geometry.y.iloc[1:],
    )

    # gapfill last period/distance values with penultimate value
    period = track.index[1:] - track.index[:-1]
    period = period.append(period[-1:])
    eye_step_distance_m = [*eye_step_distance_m, eye_step_distance_m[-1]]
    track["advection_azimuth_deg"] = [*advection_azimuth_deg, advection_azimuth_deg[-1]]

    # calculate eye speed
    track["eye_speed_ms"] = eye_step_distance_m / period.seconds.values

    # hemisphere belongs to {-1, 1}
    track["hemisphere"] = np.sign(track.geometry.y)

    adv_field: np.ndarray = np.zeros((len(track), *grid_shape), dtype=complex)
    rot_field: np.ndarray = np.zeros((len(track), *grid_shape), dtype=complex)

    for track_i, track_point in enumerate(track.itertuples()):

        adv_vector: np.complex128 = advective_vector(
            track_point.advection_azimuth_deg,
            track_point.eye_speed_ms,
            track_point.hemisphere,
        )

        adv_field[track_i, :] = np.full(grid_shape, adv_vector)

        # maximum wind speed, less advective component
        # this is the maximum tangential wind speed in the eye's non-rotating reference frame
        max_wind_speed_relative_to_eye_ms: float = track_point.max_wind_speed_ms - np.abs(adv_vector)

        rot_field[track_i, :] = rotational_field(
            longitude,  # degrees
            latitude,  # degrees
            track_point.geometry.x,  # degrees
            track_point.geometry.y,  # degrees
            track_point.radius_to_max_winds_km * 1_000,  # convert to meters
            max_wind_speed_relative_to_eye_ms,
            track_point.min_pressure_hpa * 100,  # convert to Pascals
            ENV_PRESSURE[basin] * 100,  # convert to Pascals
        )

    # sum components of wind field, (timesteps, y, x)
    wind_field: np.ndarray[complex] = adv_field + rot_field

    # take factors calculated from surface roughness of region and use to downscale speeds
    downscaled_wind_field = downscaling_factors * wind_field

    # find vector magnitude, then take max along timestep axis, giving (y, x)
    max_wind_speeds: np.ndarray[float] = np.max(np.abs(downscaled_wind_field), axis=0)

    # any dimensions with a single cell will break the plotting routines
    if 1 not in grid_shape:

        if plot_max_wind:
            plot_contours(
                max_wind_speeds,
                f"{track_id} max wind speed",
                "Wind speed [m/s]",
                os.path.join(plot_dir, f"{track_id}_max_contour.png")
            )

        if plot_animation:
            animate_track(
                downscaled_wind_field,
                track,
                os.path.join(plot_dir, f"{track_id}.gif")
            )

    return track_id, max_wind_speeds


if __name__ == "__main__":

    storm_file_path: str = snakemake.input.storm_file
    wind_grid_path: str = snakemake.input.wind_grid
    surface_roughness_path: str = snakemake.input.surface_roughness
    storm_set: set[str] = set(snakemake.params.storm_set)
    plot_max_wind: bool = snakemake.config["plot_wind"]["max_speed"]
    plot_animation: bool = snakemake.config["plot_wind"]["animation"]
    n_proc: int = snakemake.config["processes_per_parallel_job"]
    plot_dir_path: str = snakemake.output.plot_dir
    output_path: str = snakemake.output.wind_speeds
    downscale_factors_plot_path: str = snakemake.output.downscale_factors_plot

    # directory required (snakemake output)
    os.makedirs(plot_dir_path, exist_ok=True)

    logging.info("Reading tracks")
    tracks = gpd.read_parquet(storm_file_path)
    if tracks.empty:
        logging.info("No intersection between network and tracks, writing empty file.")
        cleanup(output_path, downscale_factors_plot_path)

    if storm_set:
        # if we have a storm_set, only keep the matching track_ids
        logging.info("Filtering as per storm set")
        tracks = tracks[tracks.track_id.isin(storm_set)]

    if tracks.empty:
        logging.info("No intersection between network and tracks, writing empty file.")
        cleanup(output_path, downscale_factors_plot_path)

    logging.info(f"\n{tracks}")
    grouped_tracks = tracks.groupby("track_id")

    # grid to evaluate wind speeds on, rioxarray will return midpoints of raster cells as dims
    logging.info("Reading wind evaluation grid")
    grid: xr.DataArray = rioxarray.open_rasterio(wind_grid_path)
    logging.info(f"\n{grid}")

    logging.info("Calculating downscaling factors from surface roughness raster")

    # surface roughness raster for downscaling winds with
    surface_roughness_raster: xr.DataArray = rioxarray.open_rasterio(surface_roughness_path)

    # (1, y, x) where 1 is the number of surface roughness bands
    # select the only value in the band dimension
    surface_roughness: np.ndarray = surface_roughness_raster[0].values

    # calculate factors to scale wind speeds from gradient-level to surface level,
    # taking into account surface roughness as defined by the raster
    downscaling_factors = power_law_scale_factors(
        surface_roughness,
        SURFACE_LEVEL_METRES,
        GRADIENT_LEVEL_METRES
    )

    logging.info("Mapping downscaling factors and saving to disk")
    plot_downscale_factors(
        downscaling_factors,
        "Wind downscaling factors",
        downscale_factors_plot_path
    )

    # track is a tuple of track_id and the tracks subset, we only want the latter
    args = ((track[1], grid.x, grid.y, downscaling_factors, plot_max_wind, plot_animation, plot_dir_path) for track in grouped_tracks)

    logging.info(f"Estimating wind fields for {len(grouped_tracks)} storm tracks")
    max_wind_speeds: list[str, np.ndarray] = []
    if n_proc > 1:
        with multiprocessing.Pool(processes=n_proc) as pool:
            max_wind_speeds = pool.starmap(process_track, args)
    else:
        for arg in args:
            max_wind_speeds.append(process_track(*arg))

    # sort by track_id so we have a reproducible order even after multiprocessing
    max_wind_speeds = sorted(max_wind_speeds, key=lambda pair: pair[0])

    logging.info("Saving maximum wind speeds to disk")
    track_ids, fields = zip(*max_wind_speeds)

    # write to disk as netCDF with CRS
    da = xr.DataArray(
        data=np.stack(fields),
        dims=WIND_COORDS.keys(),
        coords=(
            ("event_id", list(track_ids)),
            ("latitude", grid.y.values),
            ("longitude", grid.x.values),
        ),
        attrs=dict(
            description="Maximum estimated wind speed during event",
            units="m s-1",
        ),
        name="max_wind_speed",
    )

    # TODO: write the appropriate metadata for QGIS to read this successfully
    # as it is, the lat/long values are being ignored
    # you can of course use ncview or xarray to inspect instead...
    # spatial_ref_attrs = pyproj.CRS.from_user_input(4326).to_cf()
    # da["spatial_ref"] = ((), 0, spatial_ref_attrs)
    da = da.rio.write_crs("EPSG:4326")
    encoding = {"max_wind_speed": {"zlib": True, "complevel": 9}}
    da.to_netcdf(output_path, encoding=encoding)

    logging.info("Done estimating wind fields")
