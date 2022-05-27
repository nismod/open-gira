"""Sums the count of hits on all tranmission lines for selected region, sample, storm and aggregates reconstruction cost to level


"""

import os
import sys
from shapely.geometry import shape, LineString
import fiona
import geopandas as gpd
import json
from tqdm import tqdm
import time
from geopy import distance

try:
    output_dir = snakemake.params['output_dir']
    layer_num = snakemake.params['aggregate_level']
    reconstruction_cost = snakemake.params['reconstruction_cost']
    region_eval = snakemake.params['region_eval']
    sample_eval = snakemake.params['sample_eval']
    nh_eval = snakemake.params['nh_eval']
except:
    output_dir = 'results' #sys.argv[1]
    layer_num = 1
    reconstruction_cost = 400000
    region_eval = ['NA'] #["NA"]  # list of regions to analyse (write None if none specified)
    sample_eval = ['0'] #[0]  # list of samples of ALL regions in region_eval to analyse (write None if none specified)
    nh_eval = ['0_2005_97']  # list of storms to analyse (write None if none specified)


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




def dist_km(linestring_df):   # TODO make this function non copy paste and fix the depreciaton errors. Also fix the supression errors of processing file
    """For a linestring datafrane, returns the distance in km"""
    len_km = 0

    for ii in range(len(linestring_df)):
        if type(linestring_df.iloc[ii].geometry) == type(LineString([[1,2],[3,4]])):  # check is linestring
            line_coords = list(linestring_df.iloc[ii].geometry.coords)  # extract the coordinates of a row
            line_coords = [(x[1], x[0]) for x in line_coords]  # switch lon lat
            #print(line_coords)
            if len(line_coords) >= 2:
                for jj in range(len(line_coords)-1):
                    len_km += distance.distance(line_coords[jj], line_coords[jj+1]).km  # add the km length of that row
        else:  #multistring
            for ms in range(len(linestring_df.iloc[ii]['geometry'].geoms)):
                line_coords = list(linestring_df.iloc[ii].geometry[ms].coords)  # extract the coordinates of a row
                line_coords = [(x[1], x[0]) for x in line_coords]
                #print(line_coords)
                if len(line_coords) >= 2:
                    for jj in range(len(line_coords)-1):
                        len_km += distance.distance(line_coords[jj], line_coords[jj+1]).km  # add the km length of that row
    return len_km







boxes_county_file = os.path.join(output_dir, 'power_processed', 'world_boxes_metadata.txt')
with open(boxes_county_file, 'r') as src:
    boxes_country = json.load(src)['box_country_dict']



indiv_storm_path = os.path.join(output_dir, "power_intersection", "storm_data", "individual_storms")

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

countries_relevant = set()  # list of countries which contain relevant data (rest can be ignored)
transmission_dict = dict()
transmission_lst = []
for ii, transmission_path in tqdm(enumerate(transmission_paths), desc='Iterating transmission lines', total=len(transmission_paths)):
    transmission = gpd.read_file(transmission_path)[['link', 'geometry','box_id']]
    countries_relevant = countries_relevant.union(set().union(*[set(boxes_country[box]) for box in transmission.box_id.unique()]))  # update relevant countires
    for transmission_indiv_link in set(transmission.link):
        if transmission_indiv_link in transmission_dict.keys():
            transmission_dict[transmission_indiv_link][0] += 1

    new_dict = {transmission_indiv.link: [1, transmission_indiv.geometry] for transmission_indiv in transmission.itertuples() if transmission_indiv.link not in transmission_dict.keys()}
    transmission_dict.update(new_dict)

    print('Max value: ',max([x[0] for x in transmission_dict.values()]))

transmission_master = gpd.GeoDataFrame({'link': transmission_dict.keys(), 'count_damage':[x[0] for x in transmission_dict.values()], 'geometry':[x[1] for x in transmission_dict.values()]})


folder_agg = os.path.join(output_dir, "power_output", "statistics", "aggregate")
if not os.path.exists(folder_agg):
    os.makedirs(folder_agg)

freq_hit_path = os.path.join(folder_agg, "transmission_line_frequency_hit.gpkg")
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
            code_GIDs.append(feature["properties"]["GID_1"])
    print("creating dataframe")
    code_geoms_gpd = gpd.GeoDataFrame({"geometry": code_geoms, "code": code_GIDs})



code_geoms_gpd['len'] = [len(g.geoms) for g in code_geoms_gpd['geometry']]  # include number of polygons (higher the more computationally expensive intersect function is
code_geoms_gpd_placeholder = code_geoms_gpd.copy()

code_geoms_gpd = code_geoms_gpd.sort_values('len', ascending=True)  #TODO make false


transmission_bounds = [t.bounds for t in transmission_master.geometry]
transmission_master['lon_min'], transmission_master['lat_min'], transmission_master['lon_max'], transmission_master['lat_max'] = zip(*transmission_bounds)


map_dict = {}
print(len(transmission_master))

for geom_area in tqdm(code_geoms_gpd.itertuples(), total=len(code_geoms_gpd), desc='geom_intersect'):

    minx, miny, maxx, maxy = geom_area.geometry.bounds
    s1 = time.time()
    transmission_master_bounded = transmission_master[(transmission_master.lon_min >= minx) & (transmission_master.lon_max <= maxx) & (transmission_master.lat_min >= miny) & (transmission_master.lat_max <= maxy)]  # only check within bounding box

    if len(transmission_master_bounded) >= 1:
        #bool_list = [True if t.intersects(geom_area.geometry) else False for t in tqdm(transmission_master_bounded.geometry, desc='geom bool', total=len(transmission_master_bounded))]
        bool_list = [True if t.intersects(geom_area.geometry) else False for t in transmission_master_bounded.geometry]


        overlap_tranmission = transmission_master_bounded[bool_list]

        if len(overlap_tranmission) >= 1:
            if geom_area.code not in map_dict.keys():
                map_dict[geom_area.code] = reconstruction_cost*dist_km(overlap_tranmission)
            else:
                map_dict[geom_area.code] = map_dict[geom_area.code] + reconstruction_cost*dist_km(overlap_tranmission)

            transmission_master = transmission_master[~transmission_master.link.isin(overlap_tranmission.link)]  # remove from code_geoms_gpd_transmission_master

code_geoms_gpd['reconstruction_cost'] = code_geoms_gpd['code'].map(map_dict).fillna(0)

#code_geoms_gpd = code_geoms_gpd[code_geoms_gpd['reconstruction_cost']>0]

code_geoms_gpd.to_file(freq_hit_path.replace('frequency_hit', f'reconstruction_costs'), driver='GPKG')
print('written to file')

