"""Intersects hazard with network
Method: use the gdp flow along an edge where it is sorted for target gdp in dictionaries. Faster method.
"""

import os
import pandas as pd
import geopandas as gpd
import time
from datetime import date
from shapely.geometry import LineString
import json
from tqdm import tqdm
import sys
from damage_calculator import applythreshold
import pickle

if 'linux' not in sys.platform:  # TODO
    import os
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)


try:
    region = snakemake.params["region"]
    sample = snakemake.params["sample"]
    nh = snakemake.params["nh"]
    operationfind_ = snakemake.params["op_find"]
except:
    region = 'SP'
    sample = '0'
    nh = '0_0_1'
    operationfind_ = 'True'


if operationfind_ in [True, "True", "T", 1, "1"]:
    operationfind = True
else:
    operationfind = False

def t_op(lower, upper, targets):
    """returns number of targets between the lower and upper percentage inputs"""
    if "operation_frac" in list(targets.columns):
        return len(
            [x for x in targets["operation_frac"] if upper / 100 > x >= lower / 100]
        )
    else:
        return "N/A"


def long2short(x):  # TODO copied from functions file in process
    """Converts from long notation to short notation (more efficient to work with). x is string"""
    return x.replace('intermediate_','i').replace('_box_','b').replace('target_','t').replace('source_','s').replace('conn_','c')

def short2long(x):
    """Converts from short notation back to long notation. Order is crucial to avoid incorrect replacements. x is string"""
    return x.replace('c','conn_').replace('s','source_').replace('t','target_').replace('b','_box_').replace('i','intermediate_')




def open_collapsed(box_id):
    collapse_file = os.path.join("data", "processed", "all_boxes", box_id, f'collapsed_sources_targets_{box_id}.txt')
    with open(collapse_file, 'r') as f2:  # load the file containing midpoint collapsed nodes
        collapse_dict = json.load(f2)
    return collapse_dict

def midpoint_source(source_inp):
    """Returns the dictionary of expanded sources"""
    box_id = f"box_{source_inp.split('b')[-1]}"
    collapse_dict =  open_collapsed(box_id)
    key_source_inp = f"midpoint_{source_inp.split('_s')[1].split('b')[0]}_box_{source_inp.split('b')[-1]}"  # TODO this can be improved by renaming files in _6_
    return collapse_dict[key_source_inp]["sources"]  # retrieve sources

def midpoint_target(target_inp):
    """Returns the dictionary of expanded targets"""
    box_id = f"box_{target_inp.split('b')[-1]}"
    collapse_dict =  open_collapsed(box_id)
    key_target_inp = f"midpoint_{target_inp.split('_t')[1].split('b')[0]}_box_{target_inp.split('b')[-1]}"
    return collapse_dict[key_target_inp]["targets"]  # retrieve targets

# def midpoint_expand(source_target_dict, box_id):
#     """Takes source target dictionary with gdps and expands the midpoints according to the collapsed_sources_targets file"""
#     def collapse_file_name(box_id):
#         return os.path.join('data', 'processed', 'all_boxes', box_id, f'collapsed_sources_targets_{box_id}.txt')
#
#     for key, value in source_target_dict.items():
#         if 'midpoint' in key:
#     with open()



print(f"{nh}: loading data")
# print('loading tracks')
stormfile = os.path.join(
    "data",
    "stormtracks",
    "events",
    f"STORM_DATA_IBTRACS_{region}_1000_YEARS_{sample}.txt",
)
TC = pd.read_csv(stormfile, header=None)
TC.columns = [
    "year",
    "month",
    "number",
    "step",
    "basin",
    "lat",
    "lon",
    "pressure",
    "wind",
    "radius",
    "cat",
    "landfall",
    "dis_land",
]  # https://www.nature.com/articles/s41597-020-0381-2.pdf
TC = TC[["year", "number", "lat", "lon"]]
TC["number_hur"] = (
    str(sample)
    + "_"
    + TC["year"].astype(int).astype(str)
    + "_"
    + TC["number"].astype(int).astype(str)
)

# print("Loading wind")
windfile = os.path.join(
    "data",
    "intersection",
    "storm_data",
    "all_winds",
    region,
    f"TC_r{region}_s{sample}_n{nh}.csv",
)
if not os.path.isfile(windfile):
    raise OSError(
        f"Wind file should exist but doesnt (TC_r{region}_s{sample}_n{nh}.csv)"
    )


winds_ev_all = pd.read_csv(windfile)
assert len(winds_ev_all) != 0


stats_add_master = {}

print(f"Investigating storm {nh}")
winds_ev_filtered = applythreshold(winds_ev_all)
ID_affected = list(winds_ev_filtered["ID_point"])
box_id_affected = winds_ev_filtered["box_id"].unique()

# print("- grid")
grid_data = gpd.read_file(
    os.path.join("data", "intersection", "regions", f"{region}_unit.gpkg")
)

polys_affected = grid_data[grid_data["ID_point"].isin(ID_affected)]

# set iteration variables
routeid_damaged = set() # marks which source_sink routes already damaged -> do not double count. Form: {(source1, target1), (source2, target2), ...}
totdamage = 0  # total damage (ensuring no double counts)
targetsdamaged = {}  # if the operation value of the target is desired
edges_affected = gpd.GeoDataFrame()
targets = gpd.GeoDataFrame()

start = time.time()
for jj, box_id in enumerate(box_id_affected):
    print(f"-- Examining {jj+1}/{len(box_id_affected)} -- {box_id}")
    # print("-- network edges")
    # print("-- gdp flow")


    gdp_flow_folder = os.path.join("data", "processed", "all_boxes", box_id, 'gdp_flows')
    if (
        len(os.listdir(gdp_flow_folder)) == 0
    ):  # no gdp flow could be established (usually no sources within subgraph)
        # print("-- no gdp flow - breaking loop")
        continue

    gdp_routing_file = os.path.join("data", "processed", "all_boxes", box_id, f'edge_gdp_sorted_{box_id}.txt')
    with open(gdp_routing_file, 'r') as f1:
        gdp_routing_values = json.load(f1)





    box_edges = gpd.read_file(
        os.path.join(
            "data", "processed", "all_boxes", box_id, f"network_with_gdp_{box_id}.gpkg"
        ),
        layer="edges",
    )

    # print("-- targets")
    box_targets = gpd.read_file(
        os.path.join("data", "processed", "all_boxes", box_id, f"targets_{box_id}.gpkg")
    )

    box_edges_affected = box_edges.overlay(
        polys_affected, how="intersection"
    )  # keeps edges that are affected grid points (only a part has to be in)

    edges_to_check = [long2short(edge) for edge in box_edges_affected["link"]]

    # run through each edge to examine flor damage
    for ii, edge in tqdm(  #TODO CONSIDER JUST NOTING THE EDGES -> FORM SET -> FORM UNIQUE (SOURCE,SINK) set -> FIND in .txt CORRESPONDING DAMAGE
        enumerate(edges_to_check),
        desc=f"storm_{nh}-{box_id}: affected edges",
        total=len(box_edges_affected),
    ):

        # First load data for edge
        edge_set_file = os.path.join(gdp_flow_folder, edge+".pkl")
        route_set = set()  # set of (source, sink) tuples running through that edge
        with open(edge_set_file, 'rb') as src:
            while True:
                try:
                    route_set.add(pickle.load(src))  # to open all (when writing, the append format was used)
                except:
                    break

        dmg = 0  # damage from this edge
        dmg_frac = 1  # factor to incorporate already damaged routes

        routes_eval = route_set.difference(routeid_damaged)
        for (source, target) in routes_eval:

            gdp_damaged = gdp_routing_values[source][target]

            if 'midpoint' in source and 'midpoint' in target:  # must find the (expanded) sources and targets
                target_gdp_dict = midpoint_target(target)
                source_gdp_dict = midpoint_source(source)

                for ks, vs in source_gdp_dict.items():# TODO use itertools
                    for kt, vt in target_gdp_dict.items():
                        if (ks, kt) not in routeid_damaged:  # if the (source, target) has not already been damaged
                            dmg_frac *= vs*vt
                            routeid_damaged.update({(ks, kt)})

                            if operationfind:
                                if kt in targetsdamaged:
                                    targetsdamaged[kt] += dmg_frac*gdp_damaged
                                else:
                                    targetsdamaged[kt] = dmg_frac*gdp_damaged

            elif 'midpoint' in source:  # must find the (expanded) sources
                source_gdp_dict = midpoint_source(source)
                dict_to_assess_source = {k:v for k, v in source_gdp_dict.items() if (k, target) not in routeid_damaged}
                sum_mid_source = sum([v for v in dict_to_assess_source.values()])
                dmg_frac *= sum_mid_source  # total fraction of this (midpoint, target) combination that has not yet been damaged
                routeid_damaged.update({(source_, target) for source_ in dict_to_assess_source.keys()})

                if operationfind:
                    if target in targetsdamaged:
                        targetsdamaged[target] += sum_mid_source*gdp_damaged
                    else:
                        targetsdamaged[target] = sum_mid_source*gdp_damaged

            elif 'midpoint' in target:
                target_gdp_dict = midpoint_target(target)
                dict_to_assess_target = {k:v for k, v in target_gdp_dict.items() if (source, k) not in routeid_damaged}
                dmg_frac *= sum([v for v in dict_to_assess_target.values()])  # total fraction of this (midpoint, target) combination that has not yet been damaged

                routeid_damaged.update({(source, target_) for target_ in dict_to_assess_target.keys()})  # update which (source, target) tuples have been damaged

                if operationfind:

                    for k, v in dict_to_assess_target.items():
                        if k in targetsdamaged:
                            targetsdamaged[k] += v*gdp_damaged
                        else:
                            targetsdamaged[k] = v*gdp_damaged

            else:  # no midpoints
                routeid_damaged.update({(source, target)})  # update which (source, target) tuples have been damaged


                if operationfind:
                    if target in targetsdamaged:
                        targetsdamaged[target] += gdp_damaged
                    else:
                        targetsdamaged[target] = gdp_damaged


            gdp_damaged = gdp_routing_values[source][target]
            dmg += gdp_damaged*dmg_frac



        totdamage += dmg  # add to overall storm damage

    edges_affected = edges_affected.append(box_edges_affected)
    targets = targets.append(box_targets)

print(f"- [{nh} - Master timer: ", round((time.time() - start)/60,1), "]")

print(f"{nh}: - saving")

storm_path = os.path.join(
    "data", "intersection", "storm_data", "individual_storms", region, f"storm_{nh}"
)
if not os.path.exists(storm_path):
    os.makedirs(storm_path)

if len(edges_affected) != 0:  # to prevent writing empty dataframe
    edges_affected.to_file(
        os.path.join(
            storm_path, f"edges_affected__storm_r{region}_s{sample}_n{nh}.gpkg"
        ),
        driver="GPKG",
    )
if len(polys_affected) != 0:  # to prevent writing empty dataframe
    polys_affected.to_file(
        os.path.join(
            storm_path, f"units_affected__storm_r{region}_s{sample}_n{nh}.gpkg"
        ),
        driver="GPKG",
    )


# if other boxes (in which no direct damaged units) contain targets with a reduced operation then add to targets file
if operationfind:
    box_id_unaccounted = list(
        set(
            [
                key.split('b')[-1]
                for key in targetsdamaged.keys()
                if key.split('b')[-1] not in box_id_affected
            ]
        )
    )  # get set of boxes not in box_id_affected but have target damage
    for (
        box_id_unaccounted_indiv
    ) in box_id_unaccounted:  # for each unaccounted for box, add to targets gpkg
        targets_unaccounted_indiv = gpd.read_file(
            os.path.join(
                "data",
                "processed",
                "all_boxes",
                f"box_{box_id_unaccounted_indiv}",
                f"targets_box_{box_id_unaccounted_indiv}.gpkg",
            )
        )
        targets = targets.append(targets_unaccounted_indiv)


targetsdamaged_map = {short2long(k):v for k,v in targetsdamaged.items()}  # convert back to long
if operationfind and len(targets) != 0:
    targets["gdp_loss"] = (
        targets["id"].map(targetsdamaged_map).fillna(0)
    )  # map targets that are damaged, if not in list -> no damage (.fillna(0))
    targets["operation_frac"] = round(
        (targets["gdp"] - targets["gdp_loss"]) / targets["gdp"], 10
    )

target_cols = [
    "index",
    "area_km2",
    "population",
    "population_density_at_centroid",
    "gdp_pc",
    "gdp",
    "type",
    "box_id",
    "id",
    "geometry",
]
if len(targets) == 0:  # if empty
    targets = gpd.GeoDataFrame(columns=target_cols)
    targets.loc[0, :] = [None] * len(target_cols)

targets.to_file(
    os.path.join(storm_path, f"targets__storm_r{region}_s{sample}_n{nh}.gpkg"),
    driver="GPKG",
)


# write storm track file
if len(TC) != 0:
    print(f"- writing {nh} to storm track file")
    TC_nh = TC[TC["number_hur"] == nh]
    coords_lat = list(TC_nh["lat"])
    coords_lon = list(TC_nh["lon"])
    coords = [((coords_lon[i], coords_lat[i])) for i in range(len(coords_lat))]
    storm_track = gpd.GeoDataFrame({"geometry": [LineString(coords)]})
    storm_track.to_file(
        os.path.join(storm_path, f"storm_track_r{region}_s{sample}_n{nh}.gpkg"),
        driver="GPKG",
    )


#%% add stats
today = date.today()

if operationfind:
    targets_damaged_num = len(targetsdamaged)
else:
    frac_damage = "N/A"
    targets_damaged_num = "N/A"


stats_add = {
    "Storm ID": [nh],
    "Storm Region": [region],
    "Damages (gdp)": [totdamage],
    "targets affected": [targets_damaged_num],
    "targets operational 100%>op>=75%": [t_op(75, 100, targets)],
    "targets operational 75%>op>=50%": [t_op(50, 75, targets)],
    "targets operational 50%>op>=25%": [t_op(25, 50, targets)],
    "targets operational 25%>op>0%": [t_op(0.000001, 25, targets)],
    "targets not operational (op=0%)": [t_op(-0.000001, 0.000001, targets)],
    "sim_run_date": [today.strftime("%d/%m/%Y")],
}

damagescsvpath = os.path.join(
    storm_path,
    f"storm_r{region}_s{sample}_n{nh}.txt",
)
with open(
    damagescsvpath, "w"
) as stormfile:  # open (overwrite) file for each storm year
    json.dump(stats_add, stormfile)
