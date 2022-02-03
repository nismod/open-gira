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


# TODO: remove below lines once testing complete and solely on linux
if "linux" not in sys.platform:
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)
    region = "SP"
    nh = "0_945_1"  # number of hurr to examine
    sample = 0
    operationfind = True  # Includes the operational values of the target areas (makes about 50% slower)

else:  # linux
    region, sample, nh, operationfind_ = sys.argv[1:]
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


print("loading tracks")
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
TC['number_hur'] = str(sample)+"_"+TC['year'].astype(int).astype(str)+"_"+TC['number'].astype(int).astype(str)

print("Loading wind")
windfile = os.path.join(
    "data",
    "intersection",
    "storm_data",
    "all_winds",
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

print("- grid")
grid_data = gpd.read_file(
    os.path.join("data", "intersection", "regions", f"{region}_unit.gpkg")
)

polys_affected = grid_data[grid_data["ID_point"].isin(ID_affected)]

# set iteration variables
routeid_damaged = (
    {}
)  # marks which source_sink routes already damaged -> do not double count. Form: {source:[target,target,...], source:[target, target,...],...}
totdamage = 0  # total damage (ensuring no double counts)
if operationfind:  # if the operation value of the target is desired
    targetsdamaged = {}
edges_affected = gpd.GeoDataFrame()
targets = gpd.GeoDataFrame()

start = time.time()
for jj, box_id in enumerate(box_id_affected):
    print(f"-- Examining {jj+1}/{len(box_id_affected)} -- {box_id}")
    print("-- network edges")
    print("-- gdp flow")
    with open(
        os.path.join(
            "data", "processed", "all_boxes", box_id, f"edge_gdp_sorted_{box_id}.txt"
        ),
        "r",
    ) as sortedjson:
        sorted_gdp = json.load(sortedjson)

    if (
        len(sorted_gdp) == 0
    ):  # no gdp flow could be established (usually no sources within subgraph)
        print("-- no gdp flow - breaking loop")
        continue

    box_edges = gpd.read_file(
        os.path.join(
            "data", "processed", "all_boxes", box_id, f"network_with_gdp_{box_id}.gpkg"
        ),
        layer="edges",
    )

    print("-- targets")
    box_targets = gpd.read_file(
        os.path.join("data", "processed", "all_boxes", box_id, f"targets_{box_id}.gpkg")
    )

    box_edges_affected = box_edges.overlay(
        polys_affected, how="intersection"
    )  # keeps edges that are affected grid points (only a part has to be in)




    edges_to_check = box_edges_affected["link"]

    # run through each edge to examine flor damage
    for ii, edge in tqdm(
        enumerate(edges_to_check),
        desc=f"storm_{nh}-{box_id}: affected edges",
        total=len(box_edges_affected),
    ):
        s1 = time.time()

        # preprocess the edge
        edge_dict = sorted_gdp[
            edge
        ]  # returns a dictionary of all source_sink paths running through that edge
        edge_dict = {
            k.replace("source_", "s").replace("_box_", "b").replace("target_", "t"): v
            for k, v in edge_dict.items()
        }
        sources = [item.split("_")[0] for item in edge_dict.keys()]
        sinks = [item.split("_")[1] for item in edge_dict.keys()]

        # sum damage from undamaged source_sinks
        t1 = time.time()
        dmg = 0  # total damage per box
        for sourcesink, v in edge_dict.items():
            source = sourcesink.split("_")[0]
            if source in routeid_damaged.keys():  # key exists
                if (
                    sourcesink.split("_")[1] in routeid_damaged[source]
                ):  # sourcesink already counted
                    continue  # do not add to dmg
            dmg += v
        # print('time for t1 ', time.time()-t1)


        # find target operations
        if (
            operationfind
        ):  # if the operation value of the target is desired
            t_gdp = [
                [route_id_.split("_")[1], v]
                for route_id_, v in edge_dict.items()
                if route_id_.split("_")[1]
                not in routeid_damaged.get(route_id_.split("_")[0], [])
            ]  # list of lists containing [targetnumber, gdp] only if not already damaged
            t_dam = list(set([t[0] for t in t_gdp]))  # make unique
            for t_indiv in t_dam:  # go through unique (undamaged) targets
                totdamage_t = sum(
                    [t[1] for t in t_gdp if t[0] == t_indiv]
                )  # sum damage for each (unique, undamaged) target
                tspec_name = f"target_{t_indiv[1:t_indiv.find('b')]}_box_{t_indiv[t_indiv.find('b')+1:]}"  # target name
                if tspec_name in targetsdamaged:
                    targetsdamaged[tspec_name] += totdamage_t
                else:
                    targetsdamaged[tspec_name] = totdamage_t

        # update dictionary of damaged source_sinks
        t2 = time.time()
        sources_new = [
            source for source in sources if source not in routeid_damaged.keys()
        ]
        routeid_damaged.update(
            dict(zip(sources_new, [[]] * len(sources_new)))
        )  # add new keys
        for jj, source in enumerate(sources):
            if sinks[jj] not in routeid_damaged[source]:
                #routeid_damaged[source] += [sinks[jj]]
                routeid_damaged[source] = routeid_damaged[source].copy()+[str(sinks[jj])]
        # print('time for t2 ', time.time()-t2)


        totdamage += dmg  # add to overall storm damage

    edges_affected = edges_affected.append(box_edges_affected)
    targets = targets.append(box_targets)

print(f"- [{nh} - Master timer: ", time.time() - start, "]")

print("- saving")

storm_path = os.path.join(
    "data", "intersection", "storm_data", "individual_storms", f"storm_{nh}"
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
                key[key.find("box") :]
                for key in targetsdamaged.keys()
                if key[key.find("box") :] not in box_id_affected
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
                box_id_unaccounted_indiv,
                f"targets_{box_id_unaccounted_indiv}.gpkg",
            )
        )
        targets = targets.append(targets_unaccounted_indiv)


if operationfind and len(targets) != 0:
    targets["gdp_loss"] = targets["id"].map(targetsdamaged).fillna(0)  # map targets that are damaged, if not in list -> no damage (.fillna(0))
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
    TC_nh = TC[TC['number_hur']==nh]
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
    "data",
    "intersection",
    "storm_data",
    "individual_storms",
    f"storm_{nh}",
    f"storm_r{region}_s{sample}_n{nh}.txt",
)
with open(
    damagescsvpath, "w"
) as stormfile:  # open (overwrite) file for each storm year
    json.dump(stats_add, stormfile)
