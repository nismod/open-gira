"""For each target performs parameter analysis/gathering

The storms' .gpkg fxiles are examined, rather than the combined statistics .csv file. This is much slower but has a
wider statistical applicability. Here this feature is used to save the metrics directly to the targets in a gpkg file.
A selected value (top_select) may be input to extract to the top_select % of storms (ranked by the metric examined).
If top_select = 100 then all storms are shown.

Outputs gpkg file with metrics as target features (option to select top_select % (ranked by metric) to this file)

"""

import os
import pandas as pd
from tqdm import tqdm
import geopandas as gpd

from find_targets import find_targets, avg, sm



try:
    output_dir = snakemake.params['output_dir']
    metrics_target = snakemake.params['metrics_target']
    percentile = snakemake.params['top_select']  # percentile select (in percent). Set to 100 for all
    increased_severity_sort = snakemake.params['increased_severity_sort']
    region_eval = snakemake.params['region_eval']
    sample_eval = snakemake.params['sample_eval']
    nh_eval = snakemake.params['nh_eval']
except:  # for testing only
    output_dir = 'results'
    metrics_target = ['population_without_power', 'effective_population', 'affected_population', 'mw_loss_storm', 'f_value', 'gdp_damage']
    percentile = 100
    increased_severity_sort = True
    region_eval = None #["NA"]  # list of regions to analyse (write None if none specified)
    sample_eval = None #[0]  # list of samples of ALL regions in region_eval to analyse (write None if none specified)
    nh_eval = None  # list of storms to analyse (write None if none specified)
    #raise RuntimeError("Please use snakemake to define inputs")

import sys
if 'linux' not in sys.platform:  # TODO
    import os
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)

if sample_eval != None:
    sample_eval = [str(s) if type(s)!=str else s for s in sample_eval]
    assert all(type(x)==str for x in sample_eval)

if region_eval != None:
    assert all(type(x)==str for x in region_eval)

if nh_eval != None:
    assert all(type(x)==str for x in nh_eval)



assert 0<= percentile <= 100
assert increased_severity_sort in [True, False]
increased_severity_sort_bool = str(increased_severity_sort)[0]

## Inputs ##



## ##


stat_path = os.path.join(output_dir, 'power_output', 'statistics')
csv_path = os.path.join(stat_path, 'combined_storm_statistics.csv')
stats = pd.read_csv(csv_path)



target_paths, storm_tot = find_targets(output_dir, region_eval, sample_eval, nh_eval)
assert len(target_paths) <= storm_tot


#  update target paths to include only nh in stats with the metric nonzero


stat_path_empirical = os.path.join(stat_path, 'empirical')
if not os.path.exists(stat_path_empirical):
    os.makedirs(stat_path_empirical)


metric_sortby = {'population_without_power':'population with no power (f=0)', 'effective_population':'effective population affected', 'affected_population':'population affected', 'mw_loss_storm': 'GDP losses', 'f_value': 'GDP losses', 'gdp_damage': 'GDP losses'}  # dict: sort the stats file for this metric key by its value
assert all([metric_key in metric_sortby.keys() for metric_key in metrics_target])==True  # check all keys available
metrics_target_avg = [avg(metric) for metric in metrics_target]
metrics_target_sum = [sm(metric) for metric in metrics_target]
metric_keys = metrics_target_sum+metrics_target_avg
metric_dict = dict(zip(metric_keys, [[]]*2*len(metrics_target)))


top_select_frac = int((percentile/100)*storm_tot)  # fraction
if increased_severity_sort == True:
    text_extra = " from least to most damage (i.e. percentile definition)"
else:
    text_extra = "from most to least damage (i.e. the top 'worst' storms)"
print(f"Total {storm_tot}, stats on {len(stats)} and examining {top_select_frac} {text_extra}. Length targets is {len(target_paths)}")
storm_id_metrics = {}  # dictionary {metric1: {stormids_top_quantile_for metric1...}, metric2: {stormids_top_quantile_for metric2...}, ... }
for metric in metrics_target:
    storm_id_metrics[metric] = set(stats.sort_values(metric_sortby[metric], ascending=increased_severity_sort)['Storm ID'][:top_select_frac])  # saves a set of the top selected quantile (sorted appropriately)



metric_data_base = {}  # {target1: {metric1: val, metric2: ... , geometry: geom}, target2: {... } ...}   base info, no sum or average
metric_data = {}  # wil include sums and averages of above

# FILTER TARGET PATHS FOR STORM

for jj, target_path in tqdm(enumerate(target_paths), desc='Iterating targets', total=len(target_paths)):
    storm = os.path.basename(target_path).split('_n')[-1][:-5]  # extract storm
    #print(storm)
    targets = gpd.read_file(target_path, dtype={'population':float, 'gdp_damage':float,'mw_loss_storm':float})#[['population', 'population_density_at_centroid', 'gdp', 'id', 'f_value', 'mw_loss_storm', 'gdp_damage', 'geometry']]

    # add population definitions
    zero_pop = [t.population if float(t.f_value)==0 else 0 for t in targets.itertuples()]
    targets['population_without_power'] = zero_pop
    targets['effective_population'] = (1 - targets['f_value'].astype(float)) * targets['population'].astype(float)
    aff_pop = [t.population if float(t.f_value) < 1 else 0 for t in targets.itertuples()]
    targets['affected_population'] = aff_pop


    targets['f_value'] = 1 - targets['f_value'].astype(float)  # rescale f_rescale = 1 - f (this means that the storms with no damages contain f = 0 and so, later the average is correct. After averages are taken, it is rescaled back 1 - f_rescale
    for target_indiv in targets.itertuples():
        if target_indiv.id not in metric_data.keys():
            metric_data_new = dict(zip(metrics_target, [[]]*len(metrics_target)))  # empty (sub)dict with metrics as keys
            # metric_data_new['geometry'] = target_indiv.geometry
            # metric_data_new['country'] = target_indiv.country
            metric_data_base[target_indiv.id] = metric_data_new

            metric_data[target_indiv.id] = dict({'geometry':target_indiv.geometry, 'country': target_indiv.country})  # set up for later

        for ii, metric in enumerate(metrics_target):  # add the metric data
            if storm in storm_id_metrics[metric]:  # only if in top selected quantile (sorted by metric)

                metric_data_base[target_indiv.id][metric] = metric_data_base[target_indiv.id][metric] + [getattr(target_indiv, metric)]


for target_key in metric_data.keys():  # for each target.id
    for metric in metrics_target:  # for each metric
        metric_data[target_key][avg(metric)] = sum(metric_data_base[target_key][metric])/storm_tot  # find average
        metric_data[target_key][sm(metric)] = sum(metric_data_base[target_key][metric])  # find sum

        if 'f_value' in metric:
            metric_data[target_key][sm(metric)] = None  # sum of f_value is irrelevant
            metric_data[target_key][avg(metric)] = 1 - metric_data[target_key][avg(metric)]  # rescale back to correct

targets_combined = gpd.GeoDataFrame(metric_data).T  # include transpose due to list
t_cols = list(targets_combined.columns)
t_cols.remove('geometry')  # remove from list to be forced to float
t_cols.remove('country')  # remove from list to be forced to float
targets_combined = targets_combined.astype(dict(zip(t_cols, [float]*len(t_cols))))  # force float
print(targets_combined.describe())
folder_agg = os.path.join(output_dir, "power_output", "statistics", "aggregate")
if not os.path.exists(folder_agg):
    os.makedirs(folder_agg)
targets_combined.to_file(os.path.join(folder_agg, f"targets_geo_top{percentile}{increased_severity_sort_bool}percent.gpkg"), driver='GPKG')
print('Written to file')