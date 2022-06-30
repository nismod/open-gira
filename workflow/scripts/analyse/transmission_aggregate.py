"""Sums the count of hits on all tranmission lines for selected region, sample, storm and aggregates reconstruction cost to level


"""

import os
import sys
from shapely.geometry import shape, LineString
import fiona
import geopandas as gpd
import numpy as np
import pandas as pd
import json
from tqdm import tqdm
import time
from common_functions import find_storm_files, check_srn
from geopy import distance


try:
    output_dir = snakemake.params['output_dir']
    layer_num = snakemake.params['aggregate_level']
    region_eval = snakemake.params['region_eval']
    sample_eval = snakemake.params['sample_eval']
    nh_eval = snakemake.params['nh_eval']
    thrval = snakemake.params['central_threshold']
except:
    output_dir = 'results' #sys.argv[1]
    layer_num = 1
    region_eval = ['NA'] #["NA"]  # list of regions to analyse (write None if none specified)
    sample_eval = ['0'] #[0]  # list of samples of ALL regions in region_eval to analyse (write None if none specified)
    nh_eval = ['0_2017_66']  # list of storms to analyse (write None if none specified)
    thrval = 41
    raise RuntimeError("Please use snakemake to define inputs")

if 'linux' not in sys.platform:  # TODO
    import os
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)


def eval_dist(linestring_df):
    """"Evaluate the coordinate dataframe and returns distance in km"""

    dist_tot = 0
    for ii in range(len(linestring_df)):
        if type(linestring_df.iloc[ii].geometry) == type(LineString([[1,2],[3,4]])):  # check is linestring
            line_coords = list(linestring_df.iloc[ii].geometry.coords)  # extract the coordinates of a row
            dist_tot += eval_coordist(line_coords)

        else:  #multistring
            for ms in range(len(linestring_df.iloc[ii]['geometry'].geoms)):
                line_coords = list(linestring_df.iloc[ii].geometry[ms].coords)  # extract the coordinates of a row
                dist_tot += eval_coordist(line_coords)

    return dist_tot


def eval_dist_lst(linestring_df):
    """"Evaluate the coordinate dataframe and returns distance list in km"""

    lst = []
    for ii in range(len(linestring_df)):
        dist_tot = 0
        if type(linestring_df.iloc[ii].geometry) == type(LineString([[1,2],[3,4]])):  # check is linestring
            line_coords = list(linestring_df.iloc[ii].geometry.coords)  # extract the coordinates of a row
            dist_tot = eval_coordist(line_coords)

        else:  #multistring
            for ms in range(len(linestring_df.iloc[ii]['geometry'].geoms)):
                line_coords = list(linestring_df.iloc[ii].geometry[ms].coords)  # extract the coordinates of a row
                dist_tot = eval_coordist(line_coords)

        lst.append(dist_tot)

    return lst

def eval_coordist(coords):
    """Evaluate the coordinates and returns distance in km"""
    dist = 0
    line_coords = [(x[1], x[0]) for x in coords]
    if len(line_coords) >= 2:
        for jj in range(len(line_coords)-1):
            dist += distance.distance(line_coords[jj], line_coords[jj+1]).km  # add the km length of that row
    return dist


region_eval, sample_eval, nh_eval = check_srn(region_eval, sample_eval, nh_eval)


boxes_county_file = os.path.join(output_dir, 'power_processed', 'world_boxes_metadata.txt')
with open(boxes_county_file, 'r') as src:
    boxes_country = json.load(src)['box_country_dict']



transmission_paths, storm_tot, years_tot = find_storm_files('edges', output_dir, region_eval, sample_eval, nh_eval, thrval)

return_periods_count = np.arange(1, storm_tot+1, 1)
return_periods_ = years_tot/return_periods_count
return_periods = return_periods_[::-1]


stat_path = os.path.join(output_dir, 'power_output', 'statistics')
csv_path = os.path.join(stat_path, f'combined_storm_statistics_{thrval}.csv')
stats = pd.read_csv(csv_path, keep_default_na=False)
storms_increasing_severity = stats.sort_values('reconstruction cost', ascending=True)['Storm ID']
storm_id_metrics_weighting = dict(zip(storms_increasing_severity, return_periods))  # dictionary {storm1: return_period_of_storm1, storm2: ... }

# storm_tot = len(stats)  # TODO
# #raise RuntimeError("delete above")
# print("THIS MESSAGE SHOULD BE GONE")

folder_agg = os.path.join(output_dir, "power_output", "statistics", "aggregate")
if not os.path.exists(folder_agg):
    os.makedirs(folder_agg)

freq_hit_path = os.path.join(folder_agg, "transmission_line_frequency_hit.gpkg")
recon_path = freq_hit_path.replace('frequency_hit', 'reconstruction_costs')

if len(transmission_paths) == 0:
    print("No targets could be found. Writing dummy files (for snakemake).")
    dummy = gpd.GeoDataFrame({'geometry':[None]})
    dummy.to_file(freq_hit_path, driver='GPKG')
    dummy.to_file(recon_path, driver='GPKG')

else:

    assert len(transmission_paths) <= storm_tot


    countries_relevant = set()  # list of countries which contain relevant data (rest can be ignored)
    transmission_dict = dict()
    transmission_lst = []
    for ii, transmission_path in tqdm(enumerate(transmission_paths), desc='Iterating transmission lines', total=len(transmission_paths)):

        storm = os.path.basename(transmission_path).split('_n')[-1][:-5]  # extract storm

        transmission = gpd.read_file(transmission_path)[['link', 'geometry','box_id', 'reconstruction_cost']]
        countries_relevant = countries_relevant.union(set().union(*[set(boxes_country[box]) for box in transmission.box_id.unique()]))  # update relevant countires
        for transmission_indiv_link in set(transmission.link):
            if transmission_indiv_link in transmission_dict.keys():
                transmission_dict[transmission_indiv_link][0] += 1

        weighting_factor = 1/years_tot  # return period weighting factor for annual expected. Note this is a pretty good approximation to the trapezium rule for large n
        new_dict = {transmission_indiv.link: [1, transmission_indiv.geometry,  transmission_indiv.reconstruction_cost, weighting_factor*transmission_indiv.reconstruction_cost] for transmission_indiv in transmission.itertuples() if transmission_indiv.link not in transmission_dict.keys()}
        transmission_dict.update(new_dict)

    max_val = max([x[0] for x in transmission_dict.values()])
    print('Max value: ', max_val)

    transmission_master = gpd.GeoDataFrame({'link': transmission_dict.keys(), 'count_damage':[x[0] for x in transmission_dict.values()], 'geometry':[x[1] for x in transmission_dict.values()], 'reconstruction_cost':[x[2] for x in transmission_dict.values()], 'reconstruction_cost_annual_expected':[x[3] for x in transmission_dict.values()]})

    eval_dist_lst_output = eval_dist_lst(transmission_master)

    transmission_master['reconstruction_cost_annual_expected_per_km'] = [transmission_master['reconstruction_cost_annual_expected'].iloc[jj]/eval_dist_lst_output[jj] if eval_dist_lst_output[jj] !=0 else 0 for jj in range(len(transmission_master))]

    transmission_master.to_file(freq_hit_path, driver='GPKG')

    # Then aggregate
    print('Aggregating reconstruction costs')
    with fiona.open(
        os.path.join(output_dir, "input", "adminboundaries", f"gadm36_levels.gpkg"), "r", layer=layer_num
    ) as src_code:
        code_geoms = []
        code_GIDs = []
        for feature in src_code:
            if feature["properties"]["GID_0"] in countries_relevant:  # only include search in countries that contain targets
                code_geoms.append(shape(feature["geometry"]))
                code_GIDs.append(feature["properties"][f"GID_{layer_num}"])
        print("creating dataframe")
        code_geoms_gpd = gpd.GeoDataFrame({"geometry": code_geoms, "code": code_GIDs})



    code_geoms_gpd['len'] = [len(g.geoms) for g in code_geoms_gpd['geometry']]  # include number of polygons (higher the more computationally expensive intersect function is
    code_geoms_gpd_placeholder = code_geoms_gpd.copy()

    code_geoms_gpd = code_geoms_gpd.sort_values('len', ascending=False)

    pre_len = len(transmission_master)
    transmission_master = transmission_master[~transmission_master.geometry.isna()]  # remomve null geometries
    post_len = len(transmission_master)
    if post_len/pre_len < 0.7 and pre_len > 50:  # small check
        print('Possibly issue as 30%+ of transmission lines are being removed for having no geometry')


    transmission_bounds = [t.bounds for t in transmission_master.geometry]
    transmission_master['lon_min'], transmission_master['lat_min'], transmission_master['lon_max'], transmission_master['lat_max'] = zip(*transmission_bounds)

    map_dict = {}
    map_dict_ae = {}  # annual expected
    map_dict_ll = {}  # len_lines

    for geom_area in tqdm(code_geoms_gpd.itertuples(), total=len(code_geoms_gpd), desc='geom_intersect'):

        minx, miny, maxx, maxy = geom_area.geometry.bounds
        s1 = time.time()
        transmission_master_bounded = transmission_master[(transmission_master.lon_min >= minx) & (transmission_master.lon_max <= maxx) & (transmission_master.lat_min >= miny) & (transmission_master.lat_max <= maxy)]  # only check within bounding box



        if len(transmission_master_bounded) >= 1:
            #bool_list = [True if t.intersects(geom_area.geometry) else False for t in tqdm(transmission_master_bounded.geometry, desc='geom bool', total=len(transmission_master_bounded))]
            bool_list = [True if t.intersects(geom_area.geometry) else False for t in transmission_master_bounded.geometry]

            overlap_transmission = transmission_master_bounded[bool_list]

            if len(overlap_transmission) >= 1:

                map_dict_ll[geom_area.code] = eval_dist(overlap_transmission)  # note down degree length

                if geom_area.code not in map_dict.keys() and 'reconstruction_cost' in overlap_transmission.columns.values:
                    map_dict[geom_area.code] = sum(overlap_transmission.reconstruction_cost*overlap_transmission.count_damage)
                    map_dict_ae[geom_area.code] = sum(overlap_transmission.reconstruction_cost_annual_expected*overlap_transmission.count_damage)
                else:
                    map_dict[geom_area.code] = map_dict[geom_area.code] + sum(overlap_transmission.reconstruction_cost*overlap_transmission.count_damage)
                    map_dict_ae[geom_area.code] = map_dict_ae[geom_area.code] + sum(overlap_transmission.reconstruction_cost_annual_expected*overlap_transmission.count_damage)

                transmission_master = transmission_master[~transmission_master.link.isin(overlap_transmission.link)]  # remove from code_geoms_gpd_transmission_master

    code_geoms_gpd['reconstruction_cost_sum'] = code_geoms_gpd['code'].map(map_dict).fillna(0)
    code_geoms_gpd['reconstruction_cost_avg'] = code_geoms_gpd['reconstruction_cost_sum']/storm_tot  # average over all hitting storms
    code_geoms_gpd['reconstruction_cost_annual_expected'] = code_geoms_gpd['code'].map(map_dict_ae).fillna(0)
    code_geoms_gpd['reconstruction_cost_annual_expected_fraction_normalised'] = (code_geoms_gpd['reconstruction_cost_annual_expected']/(code_geoms_gpd['code'].map(map_dict_ll))).fillna(0)  # This metric divides by the total unit length of transmission lines in that region.

    code_geoms_gpd.to_file(recon_path, driver='GPKG')
    print('written to file')

