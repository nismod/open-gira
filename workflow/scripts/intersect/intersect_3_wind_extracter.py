"""Adapted wind speed file from J Verschuur. Processes stormtracks data and returns the wind speed at each grid location."""


import numpy as np
import pandas as pd
import sys
import geopandas as gpd
import os
import time
from tqdm import tqdm
import ast
#from pathos.multiprocessing import ProcessPool, cpu_count

try:
    region = snakemake.params["region"]
    sample = snakemake.params["sample"]
    total_parallel_processes = snakemake.params["total_parallel_processes"]
    #nh_input = snakemake.params["nh_compute"]
    all_boxes = snakemake.params["all_boxes_compute"]
except:
    pass  # cant run from console without snakemake
    # region = sys.argv[1]
    # sample = sys.argv[2]
    # #nh_input = ast.l iteral_eval(sys.argv[3])
    # nh_input = list(sys.argv[3:])


min_windlocmax = 8  # minimum wind speed value to consider significant to further save   # TODO
min_windmax = 10  # minimum wind speed value to consider significant to further save   #TODO
hurr_buffer_dist = 1500  # maximum distance to consider to strom centre

if 'linux' not in sys.platform:  # TODO
    import os
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)

    region = 'NA'
    sample = '0'
    #nh_input = str(['0_0_5'])#, '0_0_1', '0_0_2', '0_0_3', '0_0_4']
    all_boxes = [f"box_{num}" for num in [884, 955, 956, 957]]#, 1028, 1029, 1030, 1031, 1102, 1103, 1104]]  # Carribean
    total_parallel_processes = 1

#nh_input = ast.literal_eval(nh_input)  # convert string to list


def t(num, t):
    print(str(num), time.time() - t)


def haversine(lon1, lat1, lon2_lst, lat2_lst):
    lon2_arr = np.array(lon2_lst)
    lat2_arr = np.array(lat2_lst)

    # convert degrees to radians
    lon1 = np.deg2rad(lon1)
    lat1 = np.deg2rad(lat1)
    lon2_arr = np.deg2rad(lon2_arr)
    lat2_arr = np.deg2rad(lat2_arr)

    # formula
    dlon = lon2_arr - lon1
    dlat = lat2_arr - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2_arr) * np.sin(dlon / 2) ** 2
    c = 2 * np.arcsin(np.sqrt(a))
    r_e = 6371
    return c * r_e


def holland_wind_field(r, wind, pressure, pressure_env, distance, lat):
    distance = distance * 1000
    r = r * 1000
    rho = 1.10
    f = np.abs(1.45842300 * 10 ** -4 * np.sin(lat))
    e = 2.71828182846
    # p_drop = 2*wind**2
    p_drop = (pressure_env - pressure) * 100
    B = rho * e * wind ** 2 / p_drop
    Vg = (
        np.sqrt(
            ((r / distance) ** B) * (wind ** 2) * np.exp(1 - (r / distance) ** B)
            + (r ** 2) * (f ** 2) / 4
        )
        - (r * f) / 2
    )
    return Vg


start = time.time()
# nodesuse = max(1, int(cpu_count()/total_parallel_processes))  # split up fairly
# if "linux" not in sys.platform:
#     nodesuse = 7

grid_box = gpd.read_file(
    os.path.join("data", "intersection", "regions", f"{region}_unit.gpkg")
)

if len(all_boxes) != 0:
    grid_box = grid_box[grid_box['box_id'].isin(all_boxes)]  # filter for boxes that are to be examined only!

totpoints = [len(grid_box)] * len(grid_box)  # for console progress purposes
idx_pts = list(range(len(totpoints)))

sample_num = [sample] * len(grid_box)

# assert type(nh_input) == list
# nh_toprocess = []
# for nh in nh_input:
#     if not os.path.isfile(
#         os.path.join(
#             "data",
#             "intersection",
#             "storm_data",
#             "all_winds",
#             region,
#             f"TC_r{region}_s{sample}_n{nh}.csv",
#         )
#     ):
#         nh_toprocess.append(nh)
# print(f"nh to process: {nh_toprocess}")

# if len(nh_toprocess) != 0:
#     nh_toprocess_num = [nh_toprocess] * len(grid_box)

print("running wind analysis...")
#pool = ProcessPool(nodes=nodesuse)

s = time.time()






list_regions = ["NI", "SA", "NA", "EP", "SI", "SP", "WP"]
environmental_pressure = [1006.5, 1014.1, 1014.1, 1008.8, 1010.6, 1008.1, 1008.3]
environ_dict = dict(zip(list_regions, environmental_pressure))


# environ_df = pd.DataFrame(
#     {"region": list_regions, "pressure": environmental_pressure}
# )

#### load in cyclone tracks for region
stormfile = os.path.join(
    "data",
    "stormtracks",
    "events",
    "STORM_DATA_IBTRACS_" + region + "_1000_YEARS_" + sample + ".txt",
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

### add unique number to cyclone data
TC["number_hur"] = (
    sample
    + "_"
    + TC["year"].astype(int).astype(str)
    + "_"
    + TC["number"].astype(int).astype(str)
)
TC["sample"] = str(sample)

#TC = TC[TC["number_hur"].isin(nh_func)]  # filter for nhs only

##### change geometry from 0-360 to -180-180
TC['lon'] = TC['lon'].apply(lambda x: x if x <= 180 else x - 360)

# ### background environmental pressure
# pressure_env = environ_df[environ_df["region"] == str(region)]

#TC_sample = TC.copy()

grid_box = grid_box.head(200)
TC_all_lst = []
#TC_all = pd.concat([TC.copy()]*len(grid_box))  # dataframe containing all data
hurr_buffer_dist = 2000
# TODO, can we drop TC columns? (some)
distance_arr = np.array([])
for unit in tqdm(grid_box.itertuples(), desc='distances', total=len(grid_box)):  # TODO kernprof me
    distance_arr = haversine(unit.longitude, unit.latitude, TC["lon"], TC["lat"])
    TC_sample = TC.copy()
    TC_sample['distance'] = distance_arr
    TC_sample["box_id"] = unit.box_id
    TC_sample["ID_point"] = unit.ID_point
    TC_sample = TC_sample[TC_sample['distance']<=hurr_buffer_dist]
    TC_all_lst.append(TC_sample)

print("concatenating")
TC_all = pd.concat(TC_all_lst, ignore_index=True)


TC_sample = None # remove later TODO
TC_all["environ_pressure"] = environ_dict[region]


print("max winds")
### get the maximum wind and minimum distance per event
max_wind = (
    TC_all.groupby(["number_hur"])["wind"]
    .max()
    .reset_index()
    .rename(columns={"wind": "wind_max"})
)
min_distance = (
    TC_all.groupby(["number_hur"])["distance"]
    .min()
    .reset_index()
    .rename(columns={"distance": "distance_min"})
)
print("performing merges")
### merge and remove based on lowest threshold for severity hurricane
TC_all = TC_all.merge(max_wind, on="number_hur")
TC_all  = TC_all[TC_all['wind_max']>min_windmax] ### only with a maximum wind of more than  # NOTE: commented so that snakemake knows all hurricane identifiers (nh)

### merge and remove based on mimum distance set
TC_all = TC_all.merge(min_distance, on="number_hur")
#TC_all = TC_all[TC_all['distance_min']<max_distance]  # NOTE: commented so that snakemake knows all hurricane identifiers (nh)


print("Holland wind field")
TC_all["wind_location"] = holland_wind_field(
    TC_all["radius"],
    TC_all["wind"],
    TC_all["pressure"],
    TC_all["environ_pressure"],
    TC_all["distance"],
    TC_all["lat"],
)

max_wind_location = (
    TC_all.groupby(["number_hur"])["wind_location"]
    .max()
    .reset_index()
    .rename(columns={"wind_location": "wind_location_max"})
)
print("merging ")
### merge and remove based on lowest threshold for severity hurricane at location
TC_all = TC_all.merge(max_wind_location, on="number_hur")
TC_all  = TC_all[TC_all['wind_location_max']>min_windlocmax] # NOTE:  so that snakemake knows all hurricane identifiers (nh)

above20ms = (
    TC_all[TC_all["wind_location"] > 20][["wind_location", "number_hur"]]
    .groupby(["number_hur"])["wind_location"]
    .count()
    .reset_index()
    .rename(columns={"wind_location": "duration_20ms"})
)
above15ms = (
    TC_all[TC_all["wind_location"] > 15][["wind_location", "number_hur"]]
    .groupby(["number_hur"])["wind_location"]
    .count()
    .reset_index()
    .rename(columns={"wind_location": "duration_15ms"})
)

### extract the maximum wind speed at point only and associated parameters
# could add other stats here
TC_all = TC_all[
    [
        "number_hur",
        "wind_location",
        "month",
        "year",
        "sample",
        "pressure",
        "distance",
        "wind",
        "wind_max",
        "ID_point",
        "box_id",

    ]]  # TODO check if need duplicates??
# ].drop_duplicates(
#     subset=["number_hur"], keep="first"
# )  # removed: .sort_values(by = 'wind_location',ascending = False)
print("finalising")
TC_all = TC_all.merge(
    above20ms, on="number_hur", how="outer"
).replace(np.nan, 0)
TC_all = TC_all.merge(
    above15ms, on="number_hur", how="outer"
).replace(np.nan, 0)
TC_all["basin"] = region



print(f"Time for grid processing: {round((time.time()-s)/60,3)} mins")



all_winds_path = os.path.join(
    "data", "intersection", "storm_data", "all_winds", region
)
if not os.path.exists(all_winds_path):
    os.makedirs(all_winds_path)
#for nh, csv_nh in output_files.groupby("number_hur"):  #
#print(f"saving {nh}")
#nh_input.remove(nh)
p = os.path.join(all_winds_path, f"TC_r{region}_s{sample}.csv")
TC_all.to_csv(p, index=False)
print(f"Total time {round((time.time()-start)/60,3)}")
