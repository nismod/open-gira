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

        # shape=(timesteps, x, y)
        wind_field: np.ndarray = wind_field_for_track(
            interpolated_track,
            wind_grid_path,
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
            plot_field(
                wind_field_maximum,
                f"{track_id} max wind speed",
                "Wind speed [m/s]",
                os.path.join(plot_dir_path, f"{track_id}_max.png")
            )

    # TODO consider output - network joined with max speeds for all storms? could be large
#    network_with_speeds = None
#    network_with_speeds.to_parquet(output_path)


def plot_field(field: np.ndarray, title: str, colorbar_label: str, file_path: str) -> None:
    """Plot a numpy array of a 2D field and save to disk."""

    fig, ax = plt.subplots(figsize=(6,6))

    # origin lower so latitude indicies increasing northwards
    img = ax.imshow(field, vmin=0, vmax=field.max(), origin="lower")

    fig.colorbar(img, ax=ax, location="right", label=colorbar_label, shrink=0.81)
    ax.set_title(title)

    fig.savefig(file_path)
    plt.close(fig)


def animate_track(wind_field: np.ndarray, track: gpd.GeoDataFrame, file_path: str) -> None:
    """Animate a storm track and save as GIF."""

    track_name, = set(track["track_id"])
    track_length, _, _ = wind_field.shape

    fig, ax = plt.subplots(figsize=(6,6))

    # origin lower so latitude indicies increasing northwards
    img = ax.imshow(np.zeros_like(wind_field[0]), vmin=0, vmax=wind_field.max(), origin="lower")
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


def holland_wind_field(
    radius_to_max_winds_km: float,
    wind_speed_ms: float,
    pressure_hpa: float,
    pressure_env_hpa: float,
    distance_m: np.ndarray,
    lat_degrees: float,
) -> np.ndarray:
    """
    Calculate wind speed at points some distance from a cyclone track point.

    See in particular Section 3.2 in Lin and Chavas (2012).

    References
    ----------
    - Lin and Chavas (2012)
      https://www.sbafla.com/method/portals/methodology/FloodJournalArticles/Lin_Chavas_JGR12_ParametricWind.pdf
    - Holland (1980)
      https://doi.org/10.1175/1520-0493(1980)108%3C1212:AAMOTW%3E2.0.CO;2
    """

    radius_to_max_winds_m = radius_to_max_winds_km * 1000
    rho = 1.10
    f = np.abs(1.45842300e-4 * np.sin(np.radians(lat_degrees)))
    e = 2.71828182846
    delta_p = (pressure_env_hpa - pressure_hpa) * 100
    # case where (pressure_env_hpa == pressure_hpa) so p_drop is zero will raise ZeroDivisionError
    B = (
        np.power(wind_speed_ms, 2) * e * rho
        + f * wind_speed_ms * radius_to_max_winds_m * e * rho
    ) / delta_p
    Vg = (
        np.sqrt(
            # case where distance_m is zero will raise ZeroDivisionError
            (
                np.power(radius_to_max_winds_m / distance_m, B)
                * B
                * delta_p
                * np.exp(0 - (radius_to_max_winds_m / distance_m) ** B)
            )
            + (np.power(radius_to_max_winds_m, 2) * np.power(f, 2) / 4)
        )
        - (f * radius_to_max_winds_m) / 2
    )
    return Vg  # , B, delta_p, f


def wind_field_for_track(track: gpd.GeoDataFrame, wind_grid_path: str, pressure_env_hpa: float) -> np.ndarray:
    """
    Evaluate wind speed at each timestep on a grid of points given a storm track.

    Arguments:
        track (gpd.GeoDataFrame): Table of storm track information
        wind_grid_path (str): Raster file path defining a grid to evaluate on
        pressure_env_hpa (float): Background pressure for this region in hPa

    Returns:
        np.ndarray: Wind speed experienced at each grid point for each timestep
            rank=3, shape=(timesteps, x, y)
    """

    raster_grid: xr.DataArray = rioxarray.open_rasterio(wind_grid_path)
    raster_shape: tuple[int, int] = (len(raster_grid.x), len(raster_grid.y))

    grid_midpoints: list[Point] = []
    for x in raster_grid.x:
        for y in raster_grid.y:
            grid_midpoints.append(Point(x, y))
    grid = gpd.GeoSeries(grid_midpoints)

    wind_speeds: np.ndarray = np.zeros((len(track), *raster_shape))
    geod_wgs84: pyproj.Geod = pyproj.CRS("epsg:4326").get_geod()

    for track_i, track_point in enumerate(track.itertuples()):
        # distances from track to evaluation points
        # https://pyproj4.github.io/pyproj/dev/api/geod.html#pyproj.Geod.inv
        _, _, distance_m = geod_wgs84.inv(
            np.full(len(grid), track_point.geometry.x),
            np.full(len(grid), track_point.geometry.y),
            grid.geometry.x,
            grid.geometry.y,
        )

        wind_speeds[track_i]: np.ndarray = holland_wind_field(
            track_point.radius_to_max_winds_km,
            track_point.max_wind_speed_ms,
            track_point.min_pressure_hpa,
            pressure_env_hpa,
            distance_m,
            track_point.geometry.y,
        ).reshape(raster_shape, order="F")  # (x, y)

    return wind_speeds


if __name__ == "__main__":
    main()
