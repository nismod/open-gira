"""Plots the empirical distribution of storms for simple statistics for a range of threshold values

WARNING: This file is mainly for testing/plotting. In order to use this file:
The selected threshold values must first (each, separately) be computed using
the snakemake workflow. The combined_statistics.csv file must be renamed to combined_statistics_{thrval}.csv (do not
move file) where {thrval} is the threshold value in m/s (integer value). A list (thrvals) must be onput in ## inputs ##
containing the {thrval} values to be examined and plotted on the same graph.

Outputs plots of metrics vs return periods (each plot contains all selected {thrval} values)

"""

import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
from find_targets import find_targets

try:
    output_dir = snakemake.params['output_dir']
except:
    output_dir = sys.argv[1]



## Inputs ##
thrvals = [17, 25, 33, 41, 49]  # threshold values for storms
metrics = ['GDP losses', 'targets with no power (f=0)', 'population affected', 'population with no power (f=0)', 'effective population affected']

## ##

stat_path = os.path.join('results', 'power_output', 'statistics')

stat_path_empirical = os.path.join(stat_path, 'empirical', 'sensitivity')
if not os.path.exists(stat_path_empirical):
    os.makedirs(stat_path_empirical)
stat_path_empirical_data = os.path.join(stat_path, 'empirical', 'data')
if not os.path.exists(stat_path_empirical_data):
    os.makedirs(stat_path_empirical_data)

_, rp_max = find_targets(output_dir, None, None, None)  # maximum return period


for jj, thrval in enumerate(thrvals):
    print(thrval)
    csv_path = os.path.join(stat_path, f'combined_storm_statistics_{thrval}.csv')
    stats = pd.read_csv(csv_path)



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



        plt.scatter(x, y, s=2, label=f'thrval = {thrval}')
        plt.xlabel('Return Period')
        plt.ylabel(metric)
        plt.title(f"Empirical - {metric}")


        #f.set_xticks([20, 200, 500])

        if jj == len(thrvals) - 1:  # Last one

            plt.show()
            plt.grid(axis='both', which='both')
            plt.xscale("log")
            plt.legend()
            plt.savefig(os.path.join(stat_path_empirical, f'empirical_{metric}__multi_thrval.png'))




#  plt.close('all')
