"""Performs an analysis on all tranmission lines for selected region, sample, storm"""

import os
import sys
import geopandas as gpd
from tqdm import tqdm
import pandas as pd

if "linux" not in sys.platform:
    # TODO remove
    import os
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)



region_eval = None #["NA"]  # list of regions to analyse (write None if none specified)
sample_eval = None #[0]  # list of samples of ALL regions in region_eval to analyse (write None if none specified)
nh_eval = None  # list of storms to analyse (write None if none specified)



indiv_storm_path = os.path.join("results", "power_intersection", "storm_data", "individual_storms")  # TODO change output_dir

samples_add = []
if region_eval == None:
    region_eval = [os.path.basename(path) for path in os.listdir(indiv_storm_path)]
    if sample_eval == None:
        for region in region_eval:
            samples_add += [{os.path.basename(path) for path in os.listdir(os.path.join(indiv_storm_path, region))}]

        sample_eval = samples_add[0]
        for samples_test in samples_add[1:]:
            sample_eval = sample_eval.union(samples_test)  # only keep overlapping samples. Should not remove samples if correctly entered in config.

        if len(sample_eval) != len(samples_add[0]):
            print(f'Not all samples could be evaluated due to overlap. Checking samples: {sample_eval}')

transmission_paths = []
for region in region_eval:
    for sample in sample_eval:
        storms = [os.path.basename(path) for path in os.listdir(os.path.join(indiv_storm_path, region, sample))]
        transmission_paths += [os.path.join(indiv_storm_path, region, sample, file, f"edges_affected__storm_r{region}_s{sample}_n{file[6:]}.gpkg") for file in storms if os.path.isfile(os.path.join(indiv_storm_path, region, sample, file, f"edges_affected__storm_r{region}_s{sample}_n{file[6:]}.gpkg"))]  # add only if targets gpkg file exists





if len(transmission_paths) == 0:
    print('No target paths...!')

transmission_dict = dict()
transmission_lst = []
for ii, transmission_path in tqdm(enumerate(transmission_paths), desc='Iterating transmission lines', total=len(transmission_paths)):
    transmission = gpd.read_file(transmission_path)[['link', 'geometry']]
    for transmission_indiv_link in set(transmission.link):
        if transmission_indiv_link in transmission_dict.keys():
            transmission_dict[transmission_indiv_link][0]+=1

    new_dict = {transmission_indiv.link: [1, transmission_indiv.geometry] for transmission_indiv in transmission.itertuples() if transmission_indiv.link not in transmission_dict.keys()}
    transmission_dict.update(new_dict)

    print(max([x[0] for x in transmission_dict.values()]))


transmission_master = gpd.GeoDataFrame({'link': transmission_dict.keys(), 'count_damage':[x[0] for x in transmission_dict.values()], 'geometry':[x[1] for x in transmission_dict.values()]})



transmission_master.to_file("temporary_transmission_file.gpkg", driver='GPKG')