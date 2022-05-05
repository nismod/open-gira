"""Plots the empirical distribution of storms for simple statistics.

Only the combined storm .csv files are examined, rather than the .gpkg target files. This is much faster but limited to the available statistics in the .txt files.
Outputs plots of metrics vs return periods
"""

import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
from find_targets import find_targets

try:
    output_dir = snakemake.params['output_dir']
    metrics = snakemake.params['metrics']
except:
    output_dir = sys.argv[1]
    metrics = sys.argv[2:]



## Inputs ##


## ##

stat_path = os.path.join(output_dir, 'power_output', 'statistics')

stat_path_empirical = os.path.join(stat_path, 'empirical')
if not os.path.exists(stat_path_empirical):
    os.makedirs(stat_path_empirical)

csv_path = os.path.join(stat_path, 'combined_storm_statistics.csv')
stats = pd.read_csv(csv_path)

_, rp_max = find_targets(output_dir, None, None, None)  # maximum return period

for ii, metric in enumerate(metrics):
    stats_sorted = stats.sort_values(metric)  # sort for damages
    stats_sorted = stats_sorted[stats_sorted[metric]!=0]  # remove zeros
    f = plt.figure(ii)
    f.set_figwidth(10)
    f.set_figheight(8)

    y = stats_sorted[metric]

    x_count = np.arange(1, rp_max, 1)
    x = rp_max/x_count  # this is how return periods are defined!!
    x = x[:len(y)][::-1]  # correct order



    plt.scatter(x, y, s=2, color='b')
    plt.xlabel('Return Period')
    plt.ylabel(metric)
    plt.title(f"Empirical - {metric}")

    plt.grid(axis='both', which='both')
    plt.xscale("log")


    plt.show()
    plt.savefig(os.path.join(stat_path_empirical, f'empirical_{metric}.png'))

print("Plotted all empirical")
#  plt.close('all')