"""
Process storm data and return the wind speed at each grid location.
"""

import os

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import pandas as pd
import pyproj
import rioxarray
from shapely.geometry.point import Point
from tqdm import tqdm
import xarray as xr

from open_gira.direct_damages import holland_wind_model


WIND_CMAP = "turbo"
MAX_SPEED = 80  # clip wind speeds above this value when plotting
WIND_PLOT_SIZE = 8  # inches width, height


def main():
    edges_split_path = snakemake.input.edges_split  # type: ignore
    # TODO filter for relevance or replace with single file for (basin/model/sample)
    storm_file_path = snakemake.input.storm_file  # type: ignore
    wind_grid_path = snakemake.input.wind_grid
    plot_wind_fields: bool = snakemake.config["plot_wind_fields"]
    plot_dir_path = snakemake.output.plot_dir
    output_path = snakemake.output.wind_speeds  # type: ignore

    # TODO check config to restrict analysis
    # snakemake.config.specific_storm_analysis

    network = gpd.read_parquet(edges_split_path)

    # Environmental pressure values (standard estimate of background pressure away from the
    # cyclone) are taken from the AIR hurricane model, table 3 in Butke (2012).
    # Available at:
    # https://www.air-worldwide.com/publications/air-currents/2012/the-pressures-on-increased-realism-in-tropical-cyclone-wind-speeds-through-attention-to-environmental-pressure/
    environmental_pressure = {
        "NI": 1006.5,
        "SA": 1014.1,
        "NA": 1014.1,
        "EP": 1008.8,
        "SI": 1010.6,
        "SP": 1008.1,
        "WP": 1008.3,
    }

    max_wind_speeds_by_storm: dict[str, np.array] = {}
    tracks: pd.core.groupby.generic.DataFrameGroupBy = gpd.read_parquet(storm_file_path).groupby("track_id")
    for track_id, track in tqdm(tracks):

        # basin of first record for storm track (storm genesis for synthetic tracks)
        basin: str = track.iloc[0, track.columns.get_loc("basin_id")]

        interpolated_track = interpolate_track(track)

        # grid to evaluate wind speeds on, rioxarray will return midpoints of raster cells
        raster_grid: xr.DataArray = rioxarray.open_rasterio(wind_grid_path)

        # return shape = (timesteps, x, y)
        wind_field: np.ndarray = wind_field_for_track(
            interpolated_track,
            raster_grid.x,
            raster_grid.y,
            environmental_pressure[basin]
        )

        # take max along timestep dimension
        wind_field_maximum: np.ndarray = wind_field.max(axis=0)
        max_wind_speeds_by_storm[track_id] = wind_field_maximum

        # directory required (snakemake output)
        os.makedirs(plot_dir_path, exist_ok=True)
        if plot_wind_fields:
            animate_track(
                wind_field,
                interpolated_track,
                os.path.join(plot_dir_path, f"{track_id}.gif")
            )
            plot_wind_field(
                wind_field_maximum,
                f"{track_id} max wind speed",
                "Wind speed [m/s]",
                os.path.join(plot_dir_path, f"{track_id}_max.png")
            )

    # TODO consider output - network joined with max speeds for all storms? could be large
#    network_with_speeds = None
#    network_with_speeds.to_parquet(output_path)


def plot_wind_field(field: np.ndarray, title: str, colorbar_label: str, file_path: str) -> None:
    """Plot a numpy array of a 2D field wind field and save to disk."""

    # origin lower so latitude indicies increasing northwards
    origin = "lower"

    fig, ax = plt.subplots(figsize=(WIND_PLOT_SIZE, WIND_PLOT_SIZE))

    img = ax.imshow(field, vmin=0, vmax=MAX_SPEED, origin=origin, cmap=WIND_CMAP)

    levels = np.linspace(0, MAX_SPEED, int((MAX_SPEED - 0) / 10) + 1)
    contour = ax.contour(field, levels, origin=origin, colors='w')

    ax.clabel(contour, fmt='%2.1f', colors='w')

    fig.colorbar(img, ax=ax, location="right", label=colorbar_label, shrink=0.81)
    stats_str = f"min: {field.min():.2f}, max: {field.max():.2f}, $\sigma:$ {field.std():.2f}"
    ax.set_title(title + "\n" + stats_str)

    fig.savefig(file_path)
    plt.close(fig)


def animate_track(wind_field: np.ndarray, track: gpd.GeoDataFrame, file_path: str) -> None:
    """Animate a storm track and save as GIF."""

    track_name, = set(track["track_id"])
    track_length, _, _ = wind_field.shape

    fig, ax = plt.subplots(figsize=(WIND_PLOT_SIZE, WIND_PLOT_SIZE))

    # origin lower so latitude indicies increasing northwards
    img = ax.imshow(np.zeros_like(wind_field[0]), vmin=0, vmax=MAX_SPEED, origin="lower", cmap=WIND_CMAP)
    fig.colorbar(img, ax=ax, location="right", label="Wind speed [m/s]", shrink=0.81)

    def init():
        img.set_data([[]])
        return img,

    def animate(i):
        img.set_array(wind_field[i])
        ax.set_title(f"{track_name} wind speed\n{i + 1} of {track_length}")
        return img,

    anim = animation.FuncAnimation(
        fig,
        animate,
        init_func=init,
        frames=track_length,
        interval=20,
        blit=True,
    )

    anim.save(file_path, fps=10)
    plt.close(fig)


def interpolate_track(track: gpd.GeoDataFrame, substeps: int = 5) -> gpd.GeoDataFrame:
    """
    Interpolate storm track data.

    Arguments:
        track (gpd.GeoDataFrame): Storm track with at least the following
            columns: geometry, min_pressure_hpa, max_wind_speed_ms,
            radius_to_max_winds_km, timestep.
        substeps (int): For a track with n initial records, interpolate to
            n * substeps records.

    Returns:
        gpd.GeoDataFrame: Track with min_pressure_hpa, max_wind_speed_ms,
            radius_to_max_winds_km and geometry columns linearly interpolated.
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

    track["x"] = track.geometry.x
    track["y"] = track.geometry.y
    track = track.drop(columns="geometry").copy()

    dfs = [track]
    columns_to_interpolate = [
        "min_pressure_hpa",
        "max_wind_speed_ms",
        "radius_to_max_winds_km",
        "x",
        "y",
    ]
    for increment in range(1, substeps):
        substep = increment / substeps
        tmp = track.copy()
        tmp[columns_to_interpolate] = np.nan
        tmp.timestep = tmp.timestep + substep
        dfs.append(tmp)

    interpolated_track = (
        pd.concat(dfs)
        .sort_values("timestep")
        .set_index("timestep")
        .interpolate(interp_method)
    )

    # throw away repeated data after last original point
    interpolated_track = interpolated_track.iloc[:(1 - substeps), :]

    interpolated_track["geometry"] = gpd.points_from_xy(
        interpolated_track.x, interpolated_track.y, crs="EPSG:4326"
    )

    return gpd.GeoDataFrame(interpolated_track).drop(columns=["x", "y"])


def wind_field_for_track(
    track: gpd.GeoDataFrame,
    x_coords: np.ndarray,
    y_coords: np.ndarray,
    pressure_env_hpa: float
) -> np.ndarray:
    """
    Evaluate wind speed at each timestep on a grid of points given a storm track.

    Arguments:
        track (gpd.GeoDataFrame): Table of storm track information
        x_coords (np.ndarray): Longitude values to construct evaluation grid
        y_coords (np.ndarray): Latitude values to construct evaluation grid
        pressure_env_hpa (float): Background pressure for this region in hPa

    Returns:
        np.ndarray: Wind speed experienced at each grid point for each timestep
            rank=3, shape=(timesteps, x, y)
    """

    X, Y = np.meshgrid(x_coords, y_coords)

    raster_shape: tuple[int, int] = X.shape  # or Y.shape
    wind_speeds: np.ndarray = np.zeros((len(track), *raster_shape))

    geod_wgs84: pyproj.Geod = pyproj.CRS("epsg:4326").get_geod()

    for track_i, track_point in enumerate(track.itertuples()):
        # distances from track to evaluation points
        # https://pyproj4.github.io/pyproj/dev/api/geod.html#pyproj.Geod.inv
        _, _, distance_m = geod_wgs84.inv(
            np.full(len(X.ravel()), track_point.geometry.x),
            np.full(len(Y.ravel()), track_point.geometry.y),
            X.ravel(),
            Y.ravel(),
        )

        wind_speeds[track_i]: np.ndarray = holland_wind_model(
            track_point.radius_to_max_winds_km * 1_000,  # convert to meters
            track_point.max_wind_speed_ms,
            track_point.min_pressure_hpa * 100,  # convert to Pascals
            pressure_env_hpa * 100,  # convert to Pascals
            distance_m.reshape(raster_shape),  # (x, y)
            track_point.geometry.y,
        )

    return wind_speeds


if __name__ == "__main__":
    main()
