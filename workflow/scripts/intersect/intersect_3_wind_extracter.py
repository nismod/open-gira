"""
Process storm data and write out maximum wind speeds at each grid pixel for
each storm.
"""

import os

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import pandas as pd
import pyproj
import rioxarray
from tqdm import tqdm
import xarray as xr

from open_gira.direct_damages import holland_wind_model


WIND_CMAP = "turbo"
MAX_SPEED = 80  # clip wind speeds above this value when plotting
WIND_PLOT_SIZE = 9  # inches width, height

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
    wind_grid_path: str= snakemake.input.wind_grid
    plot_wind_fields: bool = snakemake.config["plot_wind_fields"]
    plot_dir_path: str= snakemake.output.plot_dir
    output_path: str = snakemake.output.wind_speeds  # type: ignore

    # TODO check config to restrict analysis
    # snakemake.config.specific_storm_analysis

    if plot_wind_fields:
        # directory required (snakemake output)
        os.makedirs(plot_dir_path, exist_ok=True)

    # grid to evaluate wind speeds on, rioxarray will return midpoints of raster cells as dims
    raster_grid: xr.DataArray = rioxarray.open_rasterio(wind_grid_path)

    max_wind_speeds_by_storm: dict[str, np.array] = {}
    tracks: pd.core.groupby.generic.DataFrameGroupBy = gpd.read_parquet(storm_file_path).groupby("track_id")

    for track_id, track in tqdm(tracks):
        if len(track) > 1:
            max_wind_speeds_by_storm[track_id] = max_wind_speed_field(track, raster_grid.x, raster_grid.y, plot_wind_fields)

    # write to disk as netCDF with CRS
    # TODO: write the appropriate metadata for QGIS to read this successfully
    # as it is, the lat/long values are being ignored
    # you can of course use ncview or xarray to inspect instead...
    da = xr.DataArray(
        data=np.stack(list(max_wind_speeds_by_storm.values())),
        dims=("event_id", "lat", "long"),
        coords=(
            ("event_id", list(max_wind_speeds_by_storm.keys())),
            ("long", raster_grid.x.values),
            ("lat", raster_grid.y.values),
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


def max_wind_speed_field(track, longitude: np.ndarray, latitude: np.ndarray, plot=False) -> np.ndarray:
    """
    Interpolate a track, reconstruct the advective and rotational vector wind
    fields, sum them and take the maximum of the wind vector magnitude across
    time.

    Args:
        track (pd.core.groupby.generic.DataFrameGroupBy): Subset of DataFrame
            describing a track. Must have a temporal index and the following
            fields: `min_pressure_hpa`, `max_wind_speed_ms`,
            `radius_to_max_winds_km`.
        longitude (np.ndarray): Longitude values to construct evaluation grid
        latitude (np.ndarray): Latitude values to construct evaluation grid

    Returns:
        np.ndarray: 2D array of maximum wind speed experienced at each grid pixel
    """

    # basin of first record for storm track (storm genesis for synthetic tracks)
    basin: str = track.iloc[0, track.columns.get_loc("basin_id")]

    # interpolate track (avoid 'doughnut effect' of wind field from infrequent eye observations)
    interpolated_track: gpd.GeoDataFrame = interpolate_track(track)

    # return shape = (timesteps, y, x) for each of advective and rotational components
    # data values are complex numbers, with i -> x, j -> y
    adv_field, rot_field = wind_field_components(
        interpolated_track,
        longitude,
        latitude,
        ENV_PRESSURE[basin]
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
            os.path.join(plot_dir_path, f"{track_id}_max_contour.png")
        )
        animate_track(
            wind_field,
            interpolated_track,
            os.path.join(plot_dir_path, f"{track_id}.gif")
        )

    return max_wind_speeds


def plot_quivers(field: np.ndarray, title: str, colorbar_label: str, file_path: str) -> None:
    """Plot a 2D numpy array of complex numbers as vector field and save to disk."""

    # origin lower so latitude indicies increasing northwards
    origin = "lower"

    fig, ax = plt.subplots(figsize=(WIND_PLOT_SIZE, WIND_PLOT_SIZE))

    ax.quiver(field.real, field.imag, angles='xy', scale=20, color='white')

    mag = np.abs(field)
    img = ax.imshow(mag, vmin=0, vmax=MAX_SPEED, origin=origin, cmap=WIND_CMAP)
    fig.colorbar(img, ax=ax, location="right", label=colorbar_label, shrink=0.81)

    stats_str = fr"min: {mag.min():.2f}, max: {mag.max():.2f}, $\sigma:$ {mag.std():.2f}"
    ax.set_title(title + "\n" + stats_str)

    fig.savefig(file_path)
    plt.close(fig)

    return


def plot_contours(field: np.ndarray, title: str, colorbar_label: str, file_path: str) -> None:
    """Plot a numpy array of a 2D field wind field and save to disk."""

    # origin lower so latitude indicies increasing northwards
    origin = "lower"

    fig, ax = plt.subplots(figsize=(WIND_PLOT_SIZE, WIND_PLOT_SIZE))

    img = ax.imshow(field, vmin=0, vmax=MAX_SPEED, origin=origin, cmap=WIND_CMAP)

    levels = np.linspace(0, MAX_SPEED, int((MAX_SPEED - 0) / 10) + 1)
    contour = ax.contour(field, levels, origin=origin, colors='w')

    ax.clabel(contour, fmt='%2.1f', colors='w')

    fig.colorbar(img, ax=ax, location="right", label=colorbar_label, shrink=0.81)
    stats_str = fr"min: {field.min():.2f}, max: {field.max():.2f}, $\sigma:$ {field.std():.2f}"
    ax.set_title(title + "\n" + stats_str)

    fig.savefig(file_path)
    plt.close(fig)

    return


def animate_track(wind_field: np.ndarray, track: gpd.GeoDataFrame, file_path: str) -> None:
    """Animate a storm track and save to disk."""

    track_name, = set(track[~track["track_id"].isna()].track_id)
    track_length, _, _ = wind_field.shape

    fig, ax = plt.subplots(figsize=(WIND_PLOT_SIZE, WIND_PLOT_SIZE))

    # origin lower so latitude indicies increasing northwards
    img = ax.imshow(np.zeros_like(np.abs(wind_field[0])), vmin=0, vmax=MAX_SPEED, origin="lower", cmap=WIND_CMAP)
    fig.colorbar(img, ax=ax, location="right", label="Wind speed [m/s]", shrink=0.81)
    quiv = ax.quiver(
        np.zeros_like(wind_field[0]),
        np.zeros_like(wind_field[0]),
        angles='xy',
        color='white'
    )

    def animate(i, image, quiver):
        image.set_array(np.abs(wind_field[i]))
        quiver.set_UVC(wind_field[i].real, wind_field[i].imag)
        ax.set_title(f"{track_name} wind vector\n{i + 1} of {track_length}")
        return img, quiv

    anim = animation.FuncAnimation(
        fig,
        animate,
        frames=track_length,
        fargs=(img, quiv),
        blit=True,
    )

    anim.save(file_path, fps=7)
    plt.close(fig)

    return


def interpolate_track(track: gpd.GeoDataFrame, frequency: str = "1H", substeps: int = 3) -> gpd.GeoDataFrame:
    """
    Interpolate storm track data.

    Arguments:
        track (gpd.GeoDataFrame): Storm track with at least the following
            columns: geometry, min_pressure_hpa, max_wind_speed_ms,
            radius_to_max_winds_km, timestep. Must have a DatetimeIndex.
        frequency (str): If given track with DatetimeIndex, interpolate to
            resolution given by this pandas frequency string
        substeps (int): If given track without DatetimeIndex, interpolate
            over track.timestep with substeps per timestep

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


def wind_field_components(
    track: gpd.GeoDataFrame,
    x_coords: np.ndarray,
    y_coords: np.ndarray,
    pressure_env_hpa: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Evaluate wind speed at each timestep on a grid of points given a storm track.

    Arguments:
        track (gpd.GeoDataFrame): Table of storm track information
        x_coords (np.ndarray): Longitude values to construct evaluation grid
        y_coords (np.ndarray): Latitude values to construct evaluation grid
        pressure_env_hpa (float): Background pressure for this region in hPa

    Returns:
        tuple(np.ndarray, np.ndarray): Advective and rotational components of
            wind speed, each field of rank=3 with shape=(timesteps, x, y)
    """

    # For reconstruction of more realistic advective wind component, see section 2 of:
    # Lin, N., and D. Chavas (2012), On hurricane parametric wind and applications
    # in storm surge modeling, J. Geophys.  Res., 117, D09120, doi:10.1029/2011JD017126
    alpha = 0.56  # fractional reduction of advective wind speed
    beta = 19.2  # degrees advective winds tend to rotate past storm track (in direction of rotation)

    X, Y = np.meshgrid(x_coords, y_coords)

    raster_shape: tuple[int, int] = X.shape  # or Y.shape
    adv: np.ndarray = np.zeros((len(track), *raster_shape), dtype=complex)
    rot: np.ndarray = np.zeros((len(track), *raster_shape), dtype=complex)

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

    for track_i, track_point in enumerate(track.itertuples()):

        if track_point.geometry.y > 0:
            hemisphere = 1  # north, hence anticlockwise rotation
        else:
            hemisphere = -1  # .. and clockwise south of the equator

        # calculate advective wind field
        # bearing of advective component (storm track heading with beta correction)
        phi_a: float = np.radians(track_point.advection_azimuth_deg - hemisphere * beta)
        mag_v_a: float = track_point.eye_speed_ms * alpha
        v_a: np.complex128 = mag_v_a * np.sin(phi_a) + mag_v_a * np.cos(phi_a) * 1j

        # TODO: should advective field eventually fall off as a function of radius?
        # otherwise for large domain with a fast moving storm, risk generating
        # powerful, potentially spurious advective winds far away from the system
        adv[track_i, :] = np.full(raster_shape, v_a)

        # forward azimuth angle and distances from grid points to track eye
        grid_to_eye_azimuth_deg, _, radius_m = geod_wgs84.inv(
            X.ravel(),
            Y.ravel(),
            np.full(len(X.ravel()), track_point.geometry.x),
            np.full(len(Y.ravel()), track_point.geometry.y),
        )

        # magnitude of rotational wind component
        mag_v_r: np.ndarray = holland_wind_model(
            track_point.radius_to_max_winds_km * 1_000,  # convert to meters
            track_point.max_wind_speed_ms - mag_v_a,  # maximum wind speed, less advective component
            track_point.min_pressure_hpa * 100,  # convert to Pascals
            pressure_env_hpa * 100,  # convert to Pascals
            radius_m.reshape(raster_shape),  # (x, y)
            track_point.geometry.y,  # latitude in degrees
        )

        # azimuth of rotational component is tangent to radius, with direction set by hemisphere
        phi_r: float = np.radians(grid_to_eye_azimuth_deg.reshape(raster_shape) + hemisphere * 90)
        rot[track_i, :] = mag_v_r * np.sin(phi_r) + mag_v_r * np.cos(phi_r) * 1j

    return adv, rot


if __name__ == "__main__":
    main()
