"""
Functions for creating animations and static plots of wind fields.
"""


import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np


WIND_CMAP = "turbo"
MAX_SPEED = 80  # clip wind speeds above this value when plotting
WIND_PLOT_SIZE = 9  # inches width, height
# origin lower so latitude indicies increasing northwards
WIND_PLOT_ORIGIN = "lower"
QUIVER_SCALE=1000


def plot_quivers(field: np.ndarray, title: str, colorbar_label: str, file_path: str) -> None:
    """Plot a 2D numpy array of complex numbers as vector field and save to disk."""

    fig, ax = plt.subplots(figsize=(WIND_PLOT_SIZE, WIND_PLOT_SIZE))

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

    fig, ax = plt.subplots(figsize=(WIND_PLOT_SIZE, WIND_PLOT_SIZE))

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

    fig, ax = plt.subplots(figsize=(WIND_PLOT_SIZE, WIND_PLOT_SIZE))

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
