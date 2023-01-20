"""
Functions for creating animations and static plots of wind fields.
"""


import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np


WIND_CMAP = "turbo"
MAX_SPEED = 80  # clip wind speeds above this value when plotting
# origin lower so latitude indicies increasing northwards
WIND_PLOT_ORIGIN = "lower"
QUIVER_SCALE = 2_500  # larger values give shorter quiver arrows


def logistic_min(x: float | np.ndarray, L: float, m: float, k: float, x_0: float) -> float | np.ndarray:
    """
    Logistic function with a minimum value, m.

    Args:
        x: Input values
        L: Maximum output value
        m: Minimum output value
        k: Steepness parameter
        x_0: Location of sigmoid centre in x

    Returns:
        Output values
    """

    return m + (L - m) / (1 + np.exp(-k * (x - x_0)))


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
    img = ax.imshow(mag, vmin=0, vmax=MAX_SPEED, origin=WIND_PLOT_ORIGIN, cmap=WIND_CMAP)
    fig.colorbar(img, ax=ax, location="right", label=colorbar_label, shrink=0.81)

    stats_str = fr"min: {mag.min():.2f}, max: {mag.max():.2f}, $\sigma:$ {mag.std():.2f}"
    ax.set_title(title + "\n" + stats_str)

    fig.savefig(file_path)
    plt.close(fig)

    return


def plot_contours(field: np.ndarray, title: str, colorbar_label: str, file_path: str) -> None:
    """Plot a numpy array of a 2D field wind field and save to disk."""

    fig, ax = plt.subplots(figsize=size_plot(*field.shape[::-1]))

    img = ax.imshow(field, vmin=0, vmax=MAX_SPEED, origin=WIND_PLOT_ORIGIN, cmap=WIND_CMAP)

    levels = np.linspace(0, MAX_SPEED, int((MAX_SPEED - 0) / 10) + 1)
    contour = ax.contour(field, levels, origin=WIND_PLOT_ORIGIN, colors='w')

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
    track_length, *grid_shape = wind_field.shape

    fig, ax = plt.subplots(figsize=size_plot(*grid_shape[::-1]))

    # origin lower so latitude indicies increasing northwards
    img = ax.imshow(np.zeros(grid_shape), vmin=0, vmax=MAX_SPEED, origin=WIND_PLOT_ORIGIN, cmap=WIND_CMAP)
    fig.colorbar(img, ax=ax, location="right", label="Wind speed [m/s]", shrink=0.81)
    quiv = ax.quiver(
        np.zeros(grid_shape),
        np.zeros(grid_shape),
        angles='xy',
        scale=QUIVER_SCALE,
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
