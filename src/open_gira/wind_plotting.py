"""
Functions for creating animations and static plots of wind fields.
"""


import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation, colors, colormaps
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import xarray as xr

from open_gira.curves import logistic_min

matplotlib.use("Agg")
plt.set_loglevel("warning")

WIND_CMAP = "magma_r"
WIND_PLOT_ORIGIN = "lower"  # origin lower so latitude indicies increasing northwards
QUIVER_SCALE = 2_500  # larger values give shorter quiver arrows

# discretise the wind field chloropleth maps
MIN = 3
MAX = 75
STEP = 3
LEVELS = int((MAX - MIN) / STEP + 1)


def size_plot(i: int, j: int) -> tuple[float, float]:
    """
    Given number of cells for each dimension, calculate appropriate figure
    width and height in inches.
    """

    L = 18
    m = 6
    k = 0.05
    x_0 = 72

    return logistic_min(i, L, m, k, x_0), logistic_min(j, L, m, k, x_0)


def plot_quivers(field: np.ndarray, title: str, colorbar_label: str, file_path: str) -> None:
    """Plot a 2D numpy array of complex numbers as vector field and save to disk."""

    fig, ax = plt.subplots(figsize=size_plot(*field.shape[::-1]))

    ax.quiver(field.real, field.imag, angles='xy', scale=QUIVER_SCALE, color='white')

    mag = np.abs(field)
    img = ax.imshow(mag, vmin=0, vmax=80, origin=WIND_PLOT_ORIGIN, cmap=WIND_CMAP)
    fig.colorbar(img, ax=ax, location="right", label=colorbar_label, shrink=0.81)

    stats_str = fr"min: {mag.min():.2f}, max: {mag.max():.2f}, $\sigma:$ {mag.std():.2f}"
    ax.set_title(title + "\n" + stats_str)

    fig.savefig(file_path)
    plt.close(fig)

    return


def plot_contours(field: np.ndarray, title: str, colorbar_label: str, file_path: str) -> None:
    """Plot a numpy array of a 2D field wind field and save to disk."""

    y, x = field.shape

    fig, ax = plt.subplots(figsize=size_plot(x, y))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    ax.axes.set_aspect('equal')

    da = xr.DataArray(field.T, coords={"i": range(x), "j": range(y)})

    cmap = colormaps[WIND_CMAP]
    cmap.set_under("w")
    xr.plot.pcolormesh(
        da,
        levels=LEVELS,
        x="i",
        y="j",
        ax=ax,
        vmin=MIN,
        vmax=MAX,
        cmap=cmap,
        cbar_ax=cax
    )

    stats_str = fr"min: {field.min():.2f}, mean: {field.mean():.2f}, max: {field.max():.2f}"
    ax.set_title(title + "\n" + stats_str)

    fig.savefig(file_path)
    plt.close(fig)

    return


def plot_downscale_factors(field: np.ndarray, title: str, file_path: str) -> None:
    """Plot a numpy array of a 2D field of [0, 1] and save to disk."""

    y, x = field.shape

    fig, ax = plt.subplots(figsize=size_plot(x, y))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    ax.axes.set_aspect('equal')

    da = xr.DataArray(field.T, coords={"i": range(x), "j": range(y)})

    xr.plot.pcolormesh(da, levels=11, x="i", y="j", ax=ax, vmin=0.5, vmax=1, cmap="inferno", cbar_ax=cax)

    stats_str = fr"min: {field.min():.2f}, mean: {field.mean():.2f}, max: {field.max():.2f}"
    ax.set_title(title + "\n" + stats_str)

    fig.savefig(file_path)
    plt.close(fig)

    return


def animate_track(wind_field: np.ndarray, track: gpd.GeoDataFrame, file_path: str) -> None:
    """Animate a storm track and save to disk."""

    track_name, = set(track[~track["track_id"].isna()].track_id)
    track_length, *grid_shape = wind_field.shape

    fig, ax = plt.subplots(figsize=size_plot(*grid_shape[::-1]))
    cmap = colormaps[WIND_CMAP]
    cmap.set_under("w")

    # origin lower so latitude indicies increasing northwards
    img = ax.imshow(
        np.zeros(grid_shape),
        norm=colors.BoundaryNorm(np.linspace(MIN, MAX, LEVELS), cmap.N),
        origin=WIND_PLOT_ORIGIN,
        cmap=cmap,
    )

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(img, ax=ax, cax=cax, label="Wind speed [m/s]", extend="min")

    def animate(i, image):
        image.set_array(np.abs(wind_field[i]))
        ax.set_title(f"{track_name} wind vector\n{i + 1} of {track_length}")
        return [img]

    anim = animation.FuncAnimation(
        fig,
        animate,
        frames=track_length,
        fargs=(img,),
        blit=True,
    )

    anim.save(file_path, fps=7)
    plt.close(fig)

    return
