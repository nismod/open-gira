"""
Plot distributions of damage fraction values by asset type
"""

import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def natural_sort(to_sort):
    return sorted(to_sort, key=lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s)])


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

    damages = pd.read_parquet(damages_path)
    hazard_cols = natural_sort([c for c in damages.columns if c.startswith("hazard-")])

    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    for asset_type in set(damages.asset_type):
        x, y = near_square_layout(len(hazard_cols))
        inches_per_subplot = 5
        fig, _ = plt.subplots(
            x,
            y,
            figsize=(1.5 * x * inches_per_subplot, y * inches_per_subplot),
            layout='constrained'
        )

        bins = np.linspace(0, 1, 20)
        for i, ax in enumerate(fig.axes):

            try:
                hazard_name = hazard_cols[i]
            except IndexError:
                # empty cells in the grid, turn off axes labels, ticks, etc.
                ax.axis('off')
                continue

            ax.hist(damages.loc[damages.asset_type == asset_type, hazard_name], bins=bins)
            ax.set_yscale('log')
            ax.set_title(hazard_name)

        fig.suptitle(f"Damage fraction distributions for {asset_type}", fontsize=16)
        fig.supylabel("Frequency", fontsize=16)
        fig.supxlabel("Damage fraction", fontsize=16)
        fig.savefig(os.path.join(plots_dir, f"{asset_type}_damage_fraction.pdf"))
