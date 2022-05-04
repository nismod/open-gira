"""Takes a gpkg file and aggregates to chosen level
Ouputs a gpkg file with metrics as target features (can use QGIS to plot)
"""

import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
import fiona
from tqdm import tqdm
import geopandas as gpd
from shapely.geometry import shape
from find_targets import avg

if "linux" not in sys.platform:
    # TODO remove
    import os
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)

output_dir = "results"  # TODO

top_select = 1
layer_num = 1



metrics = ['population', 'mw_loss_storm', 'f_value', 'gdp_damage']  # will be a target feature (column in gpkg file)
metrics_avg = [avg(metric) for metric in metrics]
metric_keys = metrics+metrics_avg


with fiona.open(
    os.path.join(output_dir, "input", "adminboundaries", f"gadm36_levels.gpkg"), "r", layer=layer_num
) as src_code:
    code_geoms = []
    code_GIDs = []
    for feature in src_code:
        code_geoms.append(shape(feature["geometry"]))
        code_GIDs.append(feature["properties"][f"GID_1"])
    print("create dataframe")
    code_geoms_gpd = gpd.GeoDataFrame({"geometry": code_geoms, "code": code_GIDs})


folder_agg = os.path.join("results", "power_output", "statistics", "aggregate")  # TODO config dir
examine_file = os.path.join(folder_agg, f"targets_geo_top{int(top_select)}percent.gpkg")


quantile_file = gpd.read_file(examine_file)


map_dict = {}
for geom_area in tqdm(code_geoms_gpd.itertuples(), total=len(code_geoms_gpd), desc='geom_intersect'):
    bool_list = [True if geom.intersects(geom_area.geometry) else False for geom in quantile_file.geometry]  # TODO check intersects

    overlap_quantile = quantile_file[bool_list]
    if len(overlap_quantile) > 0:
        if geom_area.code not in map_dict.keys():
            map_dict[geom_area.code] = dict(zip(metric_keys, [[]]*len(metric_keys)))
        for metric in metrics:
            map_dict[geom_area.code][metric] = overlap_quantile[metric].sum()
            map_dict[geom_area.code][avg(metric)] = overlap_quantile[metric].mean()


for metric in metric_keys:
    map_dict_indiv = {k: v[metric] for k, v in map_dict.items()}  # {code1: metric_value1, code2: metric_value2 ... }
    code_geoms_gpd[metric] = code_geoms_gpd['code'].map(map_dict_indiv).fillna(0)

code_geoms_gpd = code_geoms_gpd[code_geoms_gpd['gdp_damage']>0]


code_geoms_gpd.to_file(examine_file.replace('percent', f'percent_aggregated_region'))