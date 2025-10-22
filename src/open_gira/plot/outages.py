"""
Functions for drawing outage maps
"""

import os

import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from shapely.geometry import Polygon
import xarray as xr


def map_outage(
    event_id: str,
    threshold: float | int,
    exposure: xr.Dataset,
    aoi: Polygon,
    aoi_targets: gpd.GeoDataFrame,
    borders: gpd.GeoDataFrame,
    track: gpd.GeoDataFrame,
) -> plt.Figure:
    """
    Plot a target outage map for a given storm and threshold.
    """

    # preprocess data
    # extract supply factor
    df = exposure.supply_factor.sel(dict(event_id=event_id, threshold=threshold))
    df = df.to_dataframe().reset_index()[["target", "supply_factor"]]
    df = df.rename(columns={"target": "id"})

    # drop targets with NaN supply_factor
    df = df[~df.supply_factor.isna()]

    # combine target information with exposure
    data = gpd.GeoDataFrame(df.merge(aoi_targets, how="inner", on="id"))
    data.geometry = data.geometry.centroid

    # categorise supply_factor
    status_cmap = {
        "DISCONNECTED": "firebrick",
        "DEGRADED": "salmon",
        "NOMINAL": "lightgrey",
        "OVERSUPPLY": "darkorchid",
    }
    status_labels = {
        "DISCONNECTED": r"Disconnected: [$s < 0.2$]",
        "DEGRADED": r"Degraded: [$0.2 \leq s < 0.8$]",
        "NOMINAL": r"Nominal: [$0.8 \leq s < 1.2$]",
        "OVERSUPPLY": r"Oversupply: [$1.2 \leq s$]",
    }
    status_bin_edges = np.array([-np.inf, 0.20, 0.80, 1.2, np.inf])
    data["connection_status"] = pd.cut(
        data.supply_factor, bins=status_bin_edges, labels=status_cmap.keys()
    )
    data["colour"] = data.connection_status.map(status_cmap)

    # create figure that is correct aspect ratio, but no larger than 16" wide or 9" tall
    min_x, min_y, max_x, max_y = aoi.bounds
    x_span = max_x - min_x
    y_span = max_y - min_y
    aspect_ratio = y_span / x_span
    max_plot_width_in = 16
    max_plot_height_in = 9

    if max_plot_width_in * aspect_ratio < max_plot_height_in:
        # tall
        x_in = max_plot_width_in
        y_in = max_plot_width_in * aspect_ratio

    else:
        # wide
        x_in = max_plot_height_in / aspect_ratio
        y_in = max_plot_height_in

    fig, ax = plt.subplots(figsize=(x_in, y_in))

    # plot landmasses and political borders
    borders.plot(ax=ax, facecolor="none", edgecolor="grey", alpha=0.5)

    # plot supply_factor
    def population_markersize(x: np.array) -> np.array:
        """Target population -> target marker size"""
        return np.log10(x) ** 4 / 10

    ax.scatter(
        data.geometry.x,
        data.geometry.y,
        c=data.colour,
        alpha=0.3,
        marker="o",
        s=population_markersize(data.population),
    )

    pop_handles = [
        # N.B. need the sqrt around the markersize for equality between scatter markers and legend markers
        Line2D(
            [],
            [],
            color=status_cmap["NOMINAL"],
            lw=0,
            marker="o",
            markersize=np.sqrt(population_markersize(p)),
            label=f"$10^{int(np.log10(p)):d}$",
        )
        for p in np.logspace(4, 7, 7 - 4 + 1)
    ]
    pop_legend = ax.legend(
        handles=pop_handles,
        title="Node population",
        loc="lower left",
        ncol=len(pop_handles),
    )
    ax.add_artist(pop_legend)

    # reverse the cmap order, so it's from oversupply to undersupply
    cmap = dict(reversed(status_cmap.items())).items()
    status_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=colour,
            label=status_labels[status],
            markersize=8,
        )
        for status, colour in cmap
        if isinstance(status, str)
    ]
    ax.legend(
        handles=status_handles,
        ncol=1,
        title="Node supply factor, $s$",
        loc="upper right",
    )

    # plot tracks with colourbar for wind speed intensity
    track_markersize = np.exp(track.category)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.01)
    ax.plot(track.geometry.x, track.geometry.y, ls="--", color="grey", alpha=1)
    track.plot(
        column="max_wind_speed_ms",
        ax=ax,
        cax=cax,
        s=track_markersize,
        alpha=0.4,
        legend=True,
    )
    cax.set_ylabel("Wind speed $[m s^{-1}]$")

    # set window to AOI (track with a buffer)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Latitude [deg]")
    ax.grid()

    if ("name" in track.columns) and ("year" in track.columns):
        (name,) = set(track.name)
        (year,) = set(track.year)
        ax.set_title(f"{event_id}: {name}, {year:d} @ {threshold:.1f} $[m s^{{-1}}]$")
    else:
        ax.set_title(f"{event_id} @ {threshold:.1f} $[m s^{{-1}}]$")

    return fig


def animate_outage_by_threshold(
    event_id: str,
    event_dir: str,
    thresholds: list[float | int],
    exposure: xr.Dataset,
    aoi: Polygon,
    aoi_targets: gpd.GeoDataFrame,
    borders: gpd.GeoDataFrame,
    track: gpd.GeoDataFrame,
) -> None:
    """
    Plot target outage maps for a given storm and set of thresholds. Create a GIF from the frames.
    """

    if not os.path.exists(event_dir):
        os.makedirs(event_dir)

    plot_paths = []
    for threshold in thresholds:

        threshold_str = f"{threshold:.2f}".replace(".", "p")
        plot_filepath = os.path.join(event_dir, f"{threshold_str}.png")

        if not os.path.exists(plot_filepath):

            # draw map at given threshold
            fig = map_outage(
                event_id, threshold, exposure, aoi, aoi_targets, borders, track
            )
            fig.savefig(plot_filepath)

        plot_paths.append(plot_filepath)

        # animate stack of maps
        animation_filename = "outage_map_by_threshold.gif"
        os.system(
            f"convert -delay 50 {' '.join(sorted(plot_paths))} {os.path.join(event_dir, animation_filename)}"
        )
