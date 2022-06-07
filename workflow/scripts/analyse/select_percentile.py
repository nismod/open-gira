"""Selects the percentile storm as specified in the config file and copies into power_output/statistics/percentile/

"""

import os
import numpy as np
import pandas as pd
from distutils.dir_util import copy_tree
from common_functions import find_storm_files, check_srn



try:
    output_dir = snakemake.params['output_dir']
    region_eval = snakemake.params['region_eval']
    sample_eval = snakemake.params['sample_eval']
    nh_eval = snakemake.params['nh_eval']
    metrics_target = snakemake.params['metrics_target']
    thrval = snakemake.params['central_threshold']
    percentile = snakemake.params['percentile']
except:  # for testing only
    # output_dir = 'results'
    # region_eval = None #["NA"]  # list of regions to analyse (write None if none specified)
    # sample_eval = None #[0]  # list of samples of ALL regions in region_eval to analyse (write None if none specified)
    # nh_eval = None  # list of storms to analyse (write None if none specified)
    # thrval = 25
    # metrics_target = ['population_without_power', 'effective_population', 'affected_population', 'mw_loss_storm', 'f_value', 'gdp_damage']
    # percentile = 99
    raise RuntimeError("Please use snakemake to define inputs")

percentile = float(percentile)
assert 0 <= percentile <= 100
region_eval, sample_eval, nh_eval = check_srn(region_eval, sample_eval, nh_eval)

target_paths, storm_tot, years_tot = find_storm_files('targets', output_dir, region_eval, sample_eval, nh_eval, thrval)

stat_path = os.path.join(output_dir, 'power_output', 'statistics')
stat_path_percentile = os.path.join(stat_path, 'percentile')
if not os.path.exists(stat_path_percentile):
    os.makedirs(stat_path_percentile)

csv_path = os.path.join(stat_path, f'combined_storm_statistics_{thrval}.csv')
stats = pd.read_csv(csv_path, keep_default_na=False)

metric_sortby = {'population_without_power':'population with no power (f=0)', 'effective_population':'effective population affected', 'affected_population':'population affected', 'mw_loss_storm': 'GDP losses', 'f_value': 'GDP losses', 'gdp_damage': 'GDP losses'}  # dict: sort the stats file for this metric key by its value
assert all([metric_key in metric_sortby.keys() for metric_key in metrics_target])==True  # check all keys available

rp_percentile = 1/((100-percentile)/100)  # calculate the return period value


for metric in metrics_target:
    x_count = np.arange(1, years_tot+1, 1)
    x_ = years_tot/x_count
    x = x_[::-1]

    idx = len(x[x<percentile])  # find index

    storm_select = stats.sort_values(metric_sortby[metric], ascending=True)  # select percentile
    storm_region = storm_select['Storm Region'].iloc[idx]
    storm_id = storm_select['Storm ID'].iloc[idx]
    storm_sample = storm_id.split('_')[0]
    storm_id_directory = os.path.join(output_dir, 'power_intersection', 'storm_data', 'individual_storms', storm_region, str(storm_sample), f"storm_{storm_id}")

    dest_folder = os.path.join(stat_path_percentile, metric)
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)

    copy_tree(storm_id_directory, dest_folder)  # copy file of selected storm
