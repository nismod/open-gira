"""Intersects hazard with network
Method: use the gdp flow along an edge where it is sorted for target gdp in dictionaries. Faster method.
"""

import os
import pandas as pd
import geopandas as gpd
import fiona
from shapely.geometry import shape, Polygon
import ast
from datetime import date

import json
from tqdm import tqdm
import sys


# TODO: remove below lines once testing complete and solely on linux
if "linux" not in sys.platform:
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)

    code = "PHL"
    region = "WP"
    nh = "0_13"  #number of hurr to examine
    sample = 0
    operationfind = True  # Includes the operational values of the target areas (makes about 38% to 55% slower)

else:  #linux
    code, region, sample, nh, operationfind_ = sys.argv[1:]
    if operationfind_ in [True, "True", "T", 1, "1"]:
        operationfind = True



def opengrid(code):
    """Returns centroids of grid"""
    with fiona.open(os.path.join("data","intersection", f"grid_{code}.gpkg"), "r") as src:
        code_geoms = []
        for feature in src:
            code_geoms.append(shape(feature['geometry']))
        country = gpd.GeoDataFrame({'geometry':code_geoms})
    return country


def opennetworkedges():
    fname = os.path.join("data","processed", "world_network_with_gdp.gpkg")
    edges = gpd.read_file(fname, layer='edges')
    return edges


def opentargets():
    fname = os.path.join("data","processed", "world_targets.gpkg")
    targets = gpd.read_file(fname)
    return targets


def gridarea(code):
    """returns a polygon for the area """
    gdf_area = opengrid(code)

    sqr = abs(gdf_area.geometry[0].x - gdf_area.geometry[1].x)/2 # TODO not ideal, better once return period maps used
    if sqr == 0:
        raise RuntimeError("Issue with polygon creating for points")

    gdf_area = gdf_area.buffer(sqr, cap_style = 3)
    gdf_area = gpd.GeoDataFrame({"geometry":gdf_area})
    gdf_area['ID_point'] = gdf_area.index
    return gdf_area


def applythreshold(winds_df):  # TODO some threshold [for now just takes top 4]
    """Returns locations where wind is 'severe' enough to affect the power lines. Based on some condition [needs work]"""

    if len(winds_df) > 4:
        thrval = winds_df.sort_values('wind_location')['wind_location'].iloc[-4]  # take highest 4 vals
        winds_df = winds_df[winds_df['wind_location'] >= thrval]

    #winds_df = winds_df[winds_df['wind_location'] > 999]  # filter all (data testing)

    return winds_df


def t_op(lower, upper):
    """returns number of targets between the lower and upper percentage inputs"""
    return len([x for x in targets['operation_frac'] if upper/100 > x >= lower/100])


print("opening data files")
print("- grid")
grid_data = gridarea(code)

print("- network edges")
network_edges = opennetworkedges()

print("- targets")
targets = opentargets()

print("- wind")
windfile = os.path.join("data", "intersection", f'TC_c{code}_r{region}_s{sample}.csv')
winds = pd.read_csv(windfile)

print("- gdp flow")
with open(os.path.join("data","processed","edge_gdp_sorted.txt"), 'r') as sortedjson:
    sorted_gdp = json.load(sortedjson)


print("\nFiltering...")
winds_ev = winds[winds['number_hur'] == nh]  # winds for nh event

winds_ev_filtered = applythreshold(winds_ev)

ID_affected = list(winds_ev_filtered['ID_point'])
polys_affected = grid_data.iloc[ID_affected]
edges_affected = network_edges.overlay(polys_affected, how='intersection')  # keeps edges that are affected grid points (only a part has to be in)


#%%
routeid_damaged = []  # marks which source_sink routes already damaged -> do not double count
totdamage = 0  # total damage (ensuring no double counts)
ii = 0

if operationfind:  # if the operation value of the target is desired
    targetsdamaged = {}


# TODO implement more efficient grid method (snail package not working for now)
for edge in tqdm(edges_affected['link'], desc='affected edges', total=len(edges_affected)):

    edge_dict = sorted_gdp[edge]  # returns a dictionary of all source_sink paths running through that edge  (could speed up by using component dictionaries)

    dmg = [v for route_id_, v in edge_dict.items() if route_id_ not in routeid_damaged]  # count gdp damage if that source_sink hasnt already been damaged (to avoid double counting)

    if operationfind:  # if the operation value of the target is desired
        tgdp = [[route_id_[-route_id_[::-1].find('_'):], v] for route_id_, v in edge_dict.items() if route_id_ not in routeid_damaged]  # list of lists containing [targetnumber, gdp] only if not already damaged
        tdam = list(set([t[0] for t in tgdp]))
        for tspec in tdam:
            totdamage_t = sum([t[1] for t in tgdp if t[0] == tspec])
            tspec_name = f'target_{tspec}'
            if tspec_name in targetsdamaged:
                targetsdamaged[tspec_name] += totdamage_t
            else:
                targetsdamaged[tspec_name] = totdamage_t

    routeid_damaged = routeid_damaged + [link for link in edge_dict.keys() if link not in routeid_damaged]  # souce_sink link damaged -> add to list
    totdamage += sum(dmg)

if len(edges_affected) == 0:  # to prevent writing empty dataframe
    edges_affected['id'] = ['None']
if len(polys_affected) == 0:  # to prevent writing empty dataframe
    polys_affected['ID_point'] = ['None']

edges_affected.to_file(os.path.join("data","intersection", "storm_data", f"storm_{nh}", f"world_edges_affected__storm_c{code}_r{region}_s{sample}_n{nh}.gpkg"), driver='GPKG')
polys_affected.to_file(os.path.join("data","intersection", "storm_data", f"storm_{nh}", f"world_region_affected__storm_c{code}_r{region}_s{sample}_n{nh}.gpkg"), driver='GPKG')

if operationfind:
    targets['gdp_loss'] = targets['id'].map(targetsdamaged).fillna(0)
    targets['operation_frac'] = round((targets['gdp']-targets['gdp_loss'])/targets['gdp'],15)

targets.to_file(os.path.join("data","intersection","storm_data", f"storm_{nh}", f"world_targets__storm_c{code}_r{region}_s{sample}_n{nh}.gpkg"), driver='GPKG')


#%% add stats
today = date.today()

stats_add = {'Storm ID':[nh], 'Country':[code], 'Storm Region':[region], 'Damages (gdp)':[totdamage], 'fraction of total targets affected':[len(targetsdamaged)/len(targets)], 'targets affected':[len(targetsdamaged)], 'targets operational 100%>op>=75%': [t_op(75,100)], 'targets operational 75%>op>=50%': [t_op(50,75)], 'targets operational 50%>op>=25%': [t_op(25,50)], 'targets operational 25%>op>=0%': [t_op(0.0001,25)], 'targets not operational (op=0%)': [t_op(-0.0001,0.0001)], 'sim_run_date':[today.strftime("%d/%m/%Y")]}
damagescsvpath = os.path.join("data", "intersection","storm_data", f"storm_{nh}", f"storm_c{code}_r{region}_s{sample}_n{nh}.txt")
with open(damagescsvpath, 'w') as stormfile:  # open (overwrite) file for each storm
    json.dump(stats_add, stormfile)


# add to combined stats csv file, check if already exists and if already in the data (if so, remove and replace with newer)
stat_csv_path = os.path.join("data", "intersection", "combined_storm_statistics.csv")
if os.path.exists(stat_csv_path):
    df_stat = pd.read_csv(stat_csv_path)
    if nh in list(df_stat['Storm ID']):
        df_stat = df_stat.drop(df_stat[df_stat['Storm ID']==nh].index)
    df_add = pd.DataFrame(stats_add)
    df_stat = df_stat.append(df_add, ignore_index=True)
else:
    df_stat = pd.DataFrame(stats_add)
df_stat.fillna("re-run me", inplace=True)  # if any empty cells put N/A
df_stat.to_csv(stat_csv_path, index=False)