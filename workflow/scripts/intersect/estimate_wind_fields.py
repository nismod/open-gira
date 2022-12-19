"""
Process storm data and write out maximum wind speeds at each grid pixel for
each storm.
"""

import os
import multiprocessing
from typing import Optional

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
import rioxarray
import xarray as xr

from open_gira.wind import holland_wind_model, advective_vector, rotational_field
from plot_wind_fields import plot_contours, animate_track


# Environmental pressure values in hPa / mbar (standard estimate of background
# pressure away from the cyclone) are taken from the AIR hurricane model, table
# 3 in Butke (2012).  Available at:
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


def main():
    storm_file_path: str = snakemake.input.storm_file  # type: ignore
    wind_grid_path: str = snakemake.input.wind_grid
    plot_wind_fields: bool = snakemake.config["plot_wind_fields"]
    plot_dir_path: str = snakemake.output.plot_dir
    output_path: str = snakemake.output.wind_speeds  # type: ignore

    # TODO: make configurable if necessary
    parallel = True

    # TODO check config to restrict analysis
    # snakemake.config.specific_storm_analysis

    # directory required (snakemake output)
    os.makedirs(plot_dir_path, exist_ok=True)

    # grid to evaluate wind speeds on, rioxarray will return midpoints of raster cells as dims
    grid: xr.DataArray = rioxarray.open_rasterio(wind_grid_path)

    tracks = gpd.read_parquet(storm_file_path).groupby("track_id")

    # track is a tuple of track_id and the tracks subset, we only want the latter
    args = [(track[1], grid.x, grid.y, plot_wind_fields, plot_dir_path) for track in tracks]

    max_wind_speeds: list[str, np.ndarray] = []
    if parallel:
        with multiprocessing.Pool() as pool:
            max_wind_speeds = pool.starmap(process_track, args)
    else:
        for arg in args:
            max_wind_speeds.append(process_track(*arg))

    # sort by track_id so we have a reproducible order even after multiprocessing
    max_wind_speeds = sorted(max_wind_speeds, key=lambda pair: pair[0])

    track_ids, fields = zip(*max_wind_speeds)

    # write to disk as netCDF with CRS
    # TODO: write the appropriate metadata for QGIS to read this successfully
    # as it is, the lat/long values are being ignored
    # you can of course use ncview or xarray to inspect instead...
    da = xr.DataArray(
        data=np.stack(fields),
        dims=("event_id", "lat", "long"),
        coords=(
            ("event_id", list(track_ids)),
            ("long", grid.x.values),
            ("lat", grid.y.values),
        ),
        attrs=dict(
            description="Maximum estimated wind speed during event",
            units="m s-1",
        ),
        name="max_wind_speed",
    )
    da = da.rio.write_crs("epsg:4326")
    da.to_netcdf(output_path)

    return


def process_track(track, longitude: np.ndarray, latitude: np.ndarray, plot: bool, plot_dir: Optional[str]) -> tuple[str, np.ndarray]:
    """
    Interpolate a track, reconstruct the advective and rotational vector wind
    fields, sum them and take the maximum of the wind vector magnitude across
    time. Optionally plot the wind fields and save to disk.

    Args:
        track (pd.core.groupby.generic.DataFrameGroupBy): Subset of DataFrame
            describing a track. Must have a temporal index and the following
            fields: `min_pressure_hpa`, `max_wind_speed_ms`,
            `radius_to_max_winds_km`.
        longitude (np.ndarray): Longitude values to construct evaluation grid
        latitude (np.ndarray): Latitude values to construct evaluation grid
        plot (bool): Whether to plot max wind field and wind field evolution.
        plot_dir (Optional[str]): Where to save optional plots.

    Returns:
        str: Track identifier
        np.ndarray: 2D array of maximum wind speed experienced at each grid pixel
    """

    track_id, = set(track.track_id)

    # we can't calculate the advective component without at least two points
    if len(track) == 1:
        return track_id, np.zeros((len(longitude), len(latitude)))

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

    grid_shape: tuple[int, int] = (len(longitude), len(latitude))
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

    # hypotenuse of u,v wind vectors, (timesteps, y, x)
    wind_speeds: np.ndarray[float] = np.abs(wind_field)

    # find in which timestep the maximum speed values is for each raster pixel
    timestep_indicies: np.ndarray[int] = np.argmax(wind_speeds, axis=0)

    # take max along timestep axis, giving (y, x)
    # equivalently, could use timestep_indicies.choose(wind_speeds), but this limited to 32 possible timesteps
    max_wind_speeds: np.ndarray[float] = np.take_along_axis(
        wind_speeds,
        # N.B. take_along_axis indicies must be same rank as choices, so reshape from e.g. (50,50) to (1,50,50)
        timestep_indicies.reshape(1, *timestep_indicies.shape),
        axis=0
    ).squeeze()  # drop the redundant first axis

    if plot:
        plot_contours(
            max_wind_speeds,
            f"{track_id} max wind speed",
            "Wind speed [m/s]",
            os.path.join(plot_dir, f"{track_id}_max_contour.png")
        )
        animate_track(
            wind_field,
            track,
            os.path.join(plot_dir, f"{track_id}.gif")
        )

    return track_id, max_wind_speeds


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


if __name__ == "__main__":
    main()
