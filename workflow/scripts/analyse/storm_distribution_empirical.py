"""Plots the empirical distribution of storms for simple statistics.

Only the combined storm .csv files are examined, rather than the .gpkg target files. This is much faster but limited to the available statistics in the .txt files.
Outputs plots of metrics vs return periods
"""

import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
from common_functions import find_storm_files, check_srn

try:
    output_dir = snakemake.params["output_dir"]
    metrics = snakemake.params["metrics"]
    region_eval = snakemake.params["region_eval"]
    sample_eval = snakemake.params["sample_eval"]
    nh_eval = snakemake.params["nh_eval"]
    central_threshold = snakemake.params["central_threshold"]
    minimum_threshold = snakemake.params["minimum_threshold"]
    maximum_threshold = snakemake.params["maximum_threshold"]
except:  # for testing only
    output_dir = "results"
    metrics = [
        "GDP losses",
        "targets with no power (f=0)",
        "population affected",
        "population with no power (f=0)",
        "effective population affected",
    ]
    region_eval = None
    sample_eval = None
    nh_eval = None
    central_threshold = 25
    minimum_threshold = 20
    maximum_threshold = 30
    raise RuntimeError("Please use snakemake to define inputs")


def stat_file(thrval):
    """Returns pandas stat file for thrval value"""
    csv_path = os.path.join(stat_path, f"combined_storm_statistics_{thrval}.csv")
    return pd.read_csv(csv_path, keep_default_na=False)


def y_vals(df, metric):
    """Returns an array of sorted points of the metric column in the df."""
    stats_sorted = df.sort_values(metric)  # sort for damages
    stats_sorted = stats_sorted[stats_sorted[metric] != 0]  # remove zeros
    y = np.array(stats_sorted[metric])
    return y


def extme(y, len_reach):
    """Extends y with zeros at front to length len_reach"""
    if len(y) > len_reach:
        raise RuntimeError("len_reach is incorrectly specified")

    if len(y) < len_reach:
        zero_extend = np.zeros(len_reach - len(y))  # extension
        if len(y) != 0:
            y = np.concatenate((zero_extend, y))  # join
        else:
            y = zero_extend
    return y


## Inputs ##


region_eval, sample_eval, nh_eval = check_srn(region_eval, sample_eval, nh_eval)


stat_path = os.path.join(output_dir, "power_output", "statistics")

stat_path_empirical = os.path.join(stat_path, "empirical")
if not os.path.exists(stat_path_empirical):
    os.makedirs(stat_path_empirical)

stat_path_empirical_data = os.path.join(stat_path_empirical, "empirical_plotting_data")
if not os.path.exists(stat_path_empirical_data):
    os.makedirs(stat_path_empirical_data)


_, storms_tot, years_tot = find_storm_files(
    "targets", output_dir, region_eval, sample_eval, nh_eval, central_threshold
)  # maximum return period (tot number of years)


def y_extend(y1, y2, y3):
    """Extends (front) to include 0s to reach len(y_i)==max_i(y_i) for all y"""
    maxlen = max(len(y1), len(y2), len(y3))
    return extme(y1, maxlen), extme(y2, maxlen), extme(y3, maxlen)


stats_min = stat_file(minimum_threshold)
stats_max = stat_file(maximum_threshold)
stats_cen = stat_file(central_threshold)

for ii, metric in enumerate(metrics):

    f = plt.figure(ii)
    f.set_figwidth(10)
    f.set_figheight(8)

    y_min_base = y_vals(stats_min, metric)
    y_max_base = y_vals(stats_max, metric)
    y_cen_base = y_vals(stats_cen, metric)

    y_min, y_max, y_cen = y_extend(y_min_base, y_max_base, y_cen_base)

    x_count = np.arange(1, len(y_cen) + 1, 1)
    x = years_tot / x_count  # this is how return periods are defined
    x = x[::-1]  # correct order

    plt.fill_between(x, y_min, y_max, color=((138 / 256, 171 / 256, 1)), alpha=0.2)
    plt.plot(
        x,
        y_min,
        color="b",
        linestyle=":",
        label=f"Minimum wind threshold ({minimum_threshold}m/s)",
    )  # plot interpolated line
    plt.plot(
        x,
        y_max,
        color="b",
        linestyle="--",
        label=f"Maximum wind threshold ({maximum_threshold}m/s)",
    )  # plot interpolated line
    plt.plot(
        x, y_cen, color="r", label=f"Central wind threshold ({central_threshold}m/s)"
    )  # plot interpolated line
    plt.scatter(x, y_cen, s=2, color="r")  # plot data points

    plt.xlabel("Return Period")
    plt.ylabel(metric)
    plt.title(f"Empirical - {metric}")

    plt.grid(axis="both", which="both")
    plt.xscale("log")

    plt.legend()
    plt.show()
    plt.savefig(os.path.join(stat_path_empirical, f"empirical_{metric}.png"))

    # save data (in case of replot)
    data_save = pd.DataFrame({"x": x, "y_min": y_min, "y_max": y_max, "y_cen": y_cen})
    data_save.to_csv(
        os.path.join(stat_path_empirical_data, f"empirical_{metric}_plotting_data.csv"),
        index=False,
    )

print("Plotted all empirical")
#  plt.close('all')
