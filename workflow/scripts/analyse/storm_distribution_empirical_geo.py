"""Plots the empirical distribution of storms for simple statistics

The storms' .gpkg fxiles are examined, rather than the combined statistics .csv file. This is much slower but has a
wider statistical applicability. Here this feature is used to save the metrics directly to the targets in a gpkg file.
A selected value (top_select) may be input to extract to the top_select % of storms (ranked by the metric examined).
If top_select = 100 then all storms are shown.

Outputs gpkg file with metrics as target features (option to select top_select % (ranked by metric) to this file)

"""

import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
import geopandas as gpd

from find_targets import find_targets, avg

if "linux" not in sys.platform:
    # TODO remove
    import os
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)



## Inputs ##


region_eval = None #["NA"]  # list of regions to analyse (write None if none specified)
sample_eval = None #[0]  # list of samples of ALL regions in region_eval to analyse (write None if none specified)
nh_eval = None  # list of storms to analyse (write None if none specified)


top_select = 1 # top quantile select (in percent). Set to 100 for all
## ##










stat_path = os.path.join('results', 'power_output', 'statistics')
csv_path = os.path.join(stat_path, 'combined_storm_statistics.csv')
stats = pd.read_csv(csv_path)



target_paths, storm_tot = find_targets("results", region_eval, sample_eval, nh_eval)  # TODO config dir
assert len(target_paths) <= storm_tot


#  update target paths to include only nh in stats with the metric nonzero


stat_path_empirical = os.path.join(stat_path, 'empirical')
if not os.path.exists(stat_path_empirical):
    os.makedirs(stat_path_empirical)
stat_path_empirical_data = os.path.join(stat_path, 'empirical', 'data')
if not os.path.exists(stat_path_empirical_data):
    os.makedirs(stat_path_empirical_data)


metrics_ylabel = {'population': 'Population affected [people]', 'mw_loss_storm': 'Total storm power loss [mw]', 'f_value': 'Operational value [-]', 'gdp_damage': 'GDP Losses [USD]'}  # dict to precisely label y axis
metrics = ['population', 'mw_loss_storm', 'f_value', 'gdp_damage']
assert all([metric_key in metrics_ylabel.keys() for metric_key in metrics])==True  # check all keys available
metrics_avg = [avg(metric) for metric in metrics]
metric_keys = metrics+metrics_avg
metric_dict = dict(zip(metric_keys, [[]]*2*len(metrics)))


top_select_frac = int((top_select/100)*storm_tot)  # top fraction
print(f"Total {storm_tot}, stats on {len(stats)} and examining {top_select_frac}. Length targets is {len(target_paths)}")
storm_id_metrics = {}  # dictionary {metric1: {stormids_top_quantile_for metric1...}, metric2: {stormids_top_quantile_for metric2...}, ... }
for metric in metrics:
    storm_id_metrics[metric] = set(stats.sort_values('population affected', ascending=False)['Storm ID'][:top_select_frac])  # TODO metric!!!  # saves a set of the top selected quantile (sorted by quantile)  #

print(storm_id_metrics)

metric_data = {}  # {target1: {metric1: val, metric2: ... , geometry: geom}, target2: {... } ...}


# FILTER TARGET PATHS FOR STORM TODO for quantiles

for jj, target_path in tqdm(enumerate(target_paths), desc='Iterating targets', total=len(target_paths)):
    storm = os.path.basename(target_path).split('_n')[-1][:-5]  # extract storm
    #print(storm)
    targets = gpd.read_file(target_path, dtype={'population':float, 'gdp_damage':float,'mw_loss_storm':float})#[['population', 'population_density_at_centroid', 'gdp', 'id', 'f_value', 'mw_loss_storm', 'gdp_damage', 'geometry']]
    targets['f_value'] = 1 - targets['f_value'].astype(float)  # rescale f_rescale = 1 - f (this means that the storms with no damages contain f = 0 and so, later the average is correct. After averages are taken, it is rescaled back 1 - f_rescale
    for target_indiv in targets.itertuples():
        if target_indiv.id not in metric_data.keys():
            metric_data_new = dict(zip(metric_keys, [[]]*len(metric_keys)))  # empty (sub)dict with metrics as keys
            metric_data_new['geometry'] = target_indiv.geometry
            metric_data[target_indiv.id] = metric_data_new

        for ii, metric in enumerate(metrics):  # add the metric data
            if storm in storm_id_metrics[metric]:  # only if in top selected quantile (sorted by metric)

                metric_data[target_indiv.id][metric] = metric_data[target_indiv.id][metric] + [getattr(target_indiv, metric)]


for target_key in metric_data.keys():
    for metric in metrics:
        metric_data[target_key][avg(metric)] = sum(metric_data[target_key][metric])/storm_tot  # find average
        metric_data[target_key][metric] = sum(metric_data[target_key][metric])  # overwrite the lists for the sum only

        if metric == 'f_value':
            metric_data[target_key][metric] = None  # sum of f_value is irrelevant
            metric_data[target_key][avg(metric)] = 1 - metric_data[target_key][avg(metric)]  # rescale back to correct

targets_combined = gpd.GeoDataFrame(metric_data).T  # include transpose due to list
t_cols = list(targets_combined.columns)
t_cols.remove('geometry')
targets_combined = targets_combined.astype(dict(zip(t_cols, [float]*len(t_cols))))
print(targets_combined.describe())
folder_agg = os.path.join("results", "power_output", "statistics", "aggregate")  # TODO config dir
if not os.path.exists(folder_agg):
    os.makedirs(folder_agg)
targets_combined.to_file(os.path.join(folder_agg, f"targets_geo_top{int(top_select)}percent.gpkg"), driver='GPKG') #
