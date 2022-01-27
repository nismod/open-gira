"""Intersects hazard with network
Method: use the gdp flow along an edge where it is sorted for target gdp in dictionaries. Faster method.
"""

import os
import pandas as pd
import geopandas as gpd
import fiona
from shapely.geometry import shape, Polygon, Point
import ast
import time
from datetime import date
import shapely.wkt as sw
from shapely.geometry import MultiLineString, LineString
import json
from tqdm import tqdm
import sys

from damage_calculator import applythreshold

operationfind = False
# TODO: remove below lines once testing complete and solely on linux
if "linux" not in sys.platform:
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)


    region = "SP"
    year = "0"  #number of hurr to examine
    sample = 0
    operationfind = False  # Includes the operational values of the target areas (makes about 38% to 55% slower)

else:  #linux
    region, sample, year, operationfind_ = sys.argv[1:]
    if operationfind_ in [True, "True", "T", 1, "1"]:
        operationfind = True


def t_op(lower, upper, targets):
    """returns number of targets between the lower and upper percentage inputs"""
    if 'operation_frac' in list(targets.columns):
        return len([x for x in targets['operation_frac'] if upper/100 > x >= lower/100])
    else:
        return 'N/A'

print('loading tracks')
stormfile = os.path.join("data", "stormtracks", "events", f"STORM_DATA_IBTRACS_{region}_1000_YEARS_{sample}.txt")
TC = pd.read_csv(stormfile, header = None)
TC.columns = ['year', 'month', 'number', 'step','basin','lat','lon','pressure','wind','radius','cat','landfall','dis_land']  # https://www.nature.com/articles/s41597-020-0381-2.pdf
TC = TC[['year','number','lat','lon']]
TC = TC[TC['year'].astype(int)==int(year)]
TC['nh'] = str(sample)+"_"+str(year)+"_"+TC['number'].astype(int).astype(str)


print("Loading wind")
windfile = os.path.join("data", "intersection", "storm_data", "all_winds", f'TC_r{region}_s{sample}_y{year}.csv')
if not os.path.isfile(windfile):
    raise OSError(f"Wind file should exist but doesnt (TC_r{region}_s{sample}_y{year}.csv)")

try:
    winds_ev_all = pd.read_csv(windfile)

except:  # file is empty
    winds_ev_all = gpd.GeoDataFrame({'number_hur':[]})
    print(f"Storms for {year} not considered.")

stats_add_master = {}
for nh in winds_ev_all['number_hur'].unique():
    print(f"- Investigating storm {nh}")
    winds_ev = winds_ev_all[winds_ev_all['number_hur']==nh]
    winds_ev_filtered = applythreshold(winds_ev)
    ID_affected = list(winds_ev_filtered['ID_point'])
    box_id_affected = winds_ev_filtered['box_id'].unique()

    print("- grid")
    grid_data = gpd.read_file(os.path.join("data", "intersection", "regions", f"{region}_grid.gpkg"))
    #grid_data['centroid'] = [sw.loads(x) for x in grid_data['centroid']]

    polys_affected = grid_data[grid_data['ID_point'].isin(ID_affected)]

    # set iteration variables
    routeid_damaged = []  # marks which source_sink routes already damaged -> do not double count
    totdamage = 0  # total damage (ensuring no double counts)
    if operationfind:  # if the operation value of the target is desired
        targetsdamaged = {}
    edges_affected = gpd.GeoDataFrame()
    targets = gpd.GeoDataFrame()


    for jj, box_id in enumerate(box_id_affected):
        print(f"-- Examining {jj+1}/{len(box_id_affected)} -- {box_id}")
        print("-- network edges")
        box_edges = gpd.read_file(os.path.join("data","processed", "all_boxes", box_id, f"network_with_gdp_{box_id}.gpkg"), layer='edges')

        print("-- gdp flow")
        with open(os.path.join("data","processed", "all_boxes", box_id, f"edge_gdp_sorted_{box_id}.txt"), 'r') as sortedjson:
            sorted_gdp = json.load(sortedjson)

        if len(sorted_gdp) == 0:  # no gdp flow could be established (usually no sources within subgraph)
            print("-- no gdp flow - breaking loop")
            continue

        print("-- targets")
        box_targets = gpd.read_file(os.path.join("data","processed", "all_boxes", box_id, f"targets_{box_id}.gpkg"))

        # TODO below make more efficient
        box_edges_affected = box_edges.overlay(polys_affected, how='intersection')  # keeps edges that are affected grid points (only a part has to be in)

        # TODO implement more efficient grid method (snail package not working for now)
        for edge in tqdm(box_edges_affected['link'], desc='affected edges', total=len(box_edges_affected)):
            s1 = time.time()

            edge_dict = sorted_gdp[edge]  # returns a dictionary of all source_sink paths running through that edge  (could speed up by using component dictionaries)
            s2 = time.time()
            #print('1', s2-s1)
            dmg = [v for route_id_, v in edge_dict.items() if route_id_ not in routeid_damaged]  # count gdp damage if that source_sink hasnt already been damaged (to avoid double counting)
            s3 = time.time()
            #print('2', s3-s2, ",dmg len", len(dmg))
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
            s4 = time.time()
            #print('3', s4-s3, ",rd len", len(routeid_damaged))
        edges_affected = edges_affected.append(box_edges_affected)
        targets = targets.append(box_targets)

    print("- saving")

    storm_path = os.path.join("data","intersection","storm_data", "individual_storms", f"storm_{nh}")
    if not os.path.exists(storm_path):
        os.makedirs(storm_path)


    if len(edges_affected) != 0:  # to prevent writing empty dataframe
        edges_affected.to_file(os.path.join(storm_path, f"world_edges_affected__storm_r{region}_s{sample}_n{nh}.gpkg"), driver='GPKG')
    if len(polys_affected) != 0:  # to prevent writing empty dataframe
        polys_affected.to_file(os.path.join(storm_path, f"world_region_affected__storm_r{region}_s{sample}_n{nh}.gpkg"), driver='GPKG')





    if operationfind and len(targets) != 0:
        targets['gdp_loss'] = targets['id'].map(targetsdamaged).fillna(0)
        targets['operation_frac'] = round((targets['gdp']-targets['gdp_loss'])/targets['gdp'],15)

    target_cols = ['index', 'area_km2', 'population', 'population_density_at_centroid', 'gdp_pc', 'gdp', 'type', 'box_id', 'id', 'geometry']
    if len(targets) == 0:  # if empty
        targets = gpd.GeoDataFrame(columns=target_cols)
        targets.loc[0,:] = [None]*len(target_cols)

    targets.to_file(os.path.join(storm_path, f"world_targets__storm_r{region}_s{sample}_n{nh}.gpkg"), driver='GPKG')


    # write storm track file
    TC_nh = TC[TC['nh']==nh]
    if len(TC_nh) != 0:
        print(f"writing {nh} to storm track file")
        coords_lat = list(TC_nh['lat'])
        coords_lon = list(TC_nh['lon'])
        coords = [((coords_lon[i],coords_lat[i])) for i in range(len(coords_lat))]
        storm_track = gpd.GeoDataFrame({'geometry': [LineString(coords)]})
        storm_track.to_file(os.path.join(storm_path, f"storm_track_{nh}.gpkg"), driver='GPKG')



    #%% add stats
    today = date.today()

    if operationfind:
        frac_damage = len(targetsdamaged)/len(targets)
        targets_damaged_num = len(targetsdamaged)
    else:
        frac_damage = 'N/A'
        targets_damaged_num = 'N/A'


    stats_add = {'Storm ID':[nh], 'Storm Region':[region], 'Damages (gdp)':[totdamage], 'fraction of total targets affected':[frac_damage], 'targets affected':[targets_damaged_num], 'targets operational 100%>op>=75%': [t_op(75,100, targets)], 'targets operational 75%>op>=50%': [t_op(50,75, targets)], 'targets operational 50%>op>=25%': [t_op(25,50, targets)], 'targets operational 25%>op>=0%': [t_op(0.0001,25, targets)], 'targets not operational (op=0%)': [t_op(-0.0001,0.0001, targets)], 'sim_run_date':[today.strftime("%d/%m/%Y")]}
    stats_add_master[nh] = stats_add

damagescsvpath = os.path.join("data", "intersection","storm_data", "damages", f"storm_r{region}_s{sample}_y{year}.txt")
with open(damagescsvpath, 'w') as stormfile:  # open (overwrite) file for each storm year
    json.dump(stats_add_master, stormfile)


# # add to combined stats csv file, check if already exists and if already in the data (if so, remove and replace with newer)
# stat_csv_path = os.path.join("data", "intersection", "combined_storm_statistics.csv")
# if os.path.exists(stat_csv_path):
#     df_stat = pd.read_csv(stat_csv_path)
#     if nh in list(df_stat['Storm ID']):
#         df_stat = df_stat.drop(df_stat[df_stat['Storm ID']==nh].index)
#     df_add = pd.DataFrame(stats_add)
#     df_stat = df_stat.append(df_add, ignore_index=True)
# else:
#     df_stat = pd.DataFrame(stats_add)
# df_stat.fillna("re-run me", inplace=True)  # if any empty cells put N/A
# df_stat.to_csv(stat_csv_path, index=False)