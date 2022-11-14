"""
Plot distributions of damage fraction values by asset type
"""

import os
import re
import sys

import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt

from open_gira.utils import natural_sort


def closest_factors(n: int) -> tuple[tuple[int, int], int]:
    """
    Find the pair of natural numbers which factorise `n` and are closest to one another.

    Arguments:
        n (int): Number to factorise

    Returns:
        tuple[int, int]: Closest natural number factors
        int: Gap between these two factors
    """

    smallest_gap = n
    closest_factors = (1, n)

    for i in range(1, n + 1):
        for j in range(1, n + 1):
            if i * j == n:
                gap = abs(i - j)
                if gap < smallest_gap:
                    smallest_gap = gap
                    closest_factors = (i, j)

    return closest_factors, smallest_gap


def near_square_layout(n: int) -> tuple[int, int]:
    """
    Compute nearly square 2D grid shape with a number of cells equal to or narrowly in excess of `n`.

    Arguments:
        n (int): Number of items to accomodate

    Returns:
        tuple[int, int]: Suggested x, y layout
    """

    MAX_EMPTY_CELLS = 3
    optimum = (1, n)

    for i in range(n, n + MAX_EMPTY_CELLS + 1):
        factors, gap = closest_factors(i)
        if gap < abs(optimum[0] - optimum[1]):
            optimum = factors

    # return larger then smaller (prefer landscape for x,y addressing)
    return tuple(sorted(optimum, reverse=True))


if __name__ == "__main__":

    try:
        damages_path = snakemake.input["damages"]
        plots_dir = snakemake.output["plots"]
    except NameError:
        damages_path, plots_dir = sys.argv[1:]

    damages = gpd.read_parquet(damages_path)

    # get length information for each edge
    if not "length_km" in damages.columns:
        damages = damages.set_crs(epsg=4326)
        damages["length_km"] = damages.to_crs(damages.estimate_utm_crs()).geometry.length / 1_000

    hazard_cols = natural_sort([c for c in damages.columns if c.startswith("hazard-")])

    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    for asset_type in set(damages.asset_type):

        # subset
        asset_damages = damages[damages.asset_type == asset_type]

        x, y = near_square_layout(len(hazard_cols))
        inches_per_subplot = 5
        fig, _ = plt.subplots(
            x,
            y,
            figsize=(1.5 * x * inches_per_subplot, y * inches_per_subplot)
        )

        n_bins = 10
        bins = np.linspace(0, 1, n_bins + 1)
        bin_width = (bins[1:] - bins[:-1]).mean()

        for i, ax in enumerate(fig.axes):

            try:
                hazard_name = hazard_cols[i]
            except IndexError:
                # empty cells in the grid, turn off axes labels, ticks, etc.
                ax.axis('off')
                continue

            # find the index of the bin each damage fraction falls into
            bin_indicies_by_row = np.digitize(asset_damages[hazard_name], bins) - 1

            # sum the lengths of binned edges
            length_by_bin = []
            for i in range(n_bins):
                segment_lengths = asset_damages.iloc[
                    np.where(bin_indicies_by_row == i)[0],
                    asset_damages.columns.get_loc("length_km")
                ]
                length_by_bin.append(segment_lengths.sum())

            # make a histogram of sorts with the length of exposed edges
            ax.bar(bins[:-1], length_by_bin, width=bin_width, align='edge')
            ax.set_yscale("log")
            ax.set_title(hazard_name)
            ax.grid()

        total_asset_length = asset_damages.length_km.sum()
        fig.suptitle(f"Damage fraction distributions for {total_asset_length:.0f} km of {asset_type}", fontsize=16)
        fig.supylabel("Length of asset at damage fraction (km)", fontsize=16)
        fig.supxlabel(f"Damage fraction", fontsize=16)
        fig.savefig(os.path.join(plots_dir, f"{asset_type}_damage_fraction.pdf"))
