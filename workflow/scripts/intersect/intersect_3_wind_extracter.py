"""Process storm data and return the wind speed at each grid location.

Adapted wind speed file from J Verschuur.

TODO: Investigate speeding up this script. N.B. ~80% of execution time is file I/O.

"""
import os
import time

import geopandas as gpd
import numpy as np
import pandas as pd
from tqdm import tqdm

try:
    region = snakemake.params["region"]  # type: ignore
    sample = snakemake.params["sample"]  # type: ignore
    all_boxes = snakemake.params["all_boxes_compute"]  # type: ignore
    nh_split = int(
        snakemake.params["memory_storm_split"]  # type: ignore
    )  # number of nh to run each iteration'
    wind_rerun = snakemake.params["wind_rerun"]  # type: ignore
    output_dir = snakemake.params["output_dir"]  # type: ignore
    stormfile = snakemake.params["storm_file"]  # type: ignore
    central_threshold = snakemake.params["central_threshold"]  # type: ignore
    minimum_threshold = snakemake.params["minimum_threshold"]  # type: ignore
    maximum_threshold = snakemake.params["maximum_threshold"]  # type: ignore
except:
    # raise RuntimeError("Snakemake parameters not found")
    region = "NA"
    sample = 0
    all_boxes = [
        f"box_{num}" for num in [884, 955, 956, 957, 1028, 1029, 1030, 1031, 1103, 1104]
    ]
    nh_split = 2500
    wind_rerun = False
    output_dir = "results"
    stormfile = "results/input/storm-ibtracs/events/constant/NA/STORM_DATA_CMCC-CM2-VHR4_NA_1000_YEARS_0_IBTRACSDELTA.txt"
    central_threshold = 43
    minimum_threshold = 39
    maximum_threshold = 47

min_windlocmax = (
    float(minimum_threshold)
    - 10  # minimum wind speed value (at unit) to consider significant to further save
)
min_windmax = (
    float(minimum_threshold) - 8
)  # minimum wind speed value (over entire storm at cyclone centre) to consider significant to further save
hurr_buffer_dist = 1300  # maximum distance to consider to storm centre


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
    lat = lat * np.pi / 180
    distance = distance * 1000
    r = r * 1000
    rho = 1.10
    f = np.abs(1.45842300 * 10**-4 * np.sin(lat))
    e = 2.71828182846
    # p_drop = 2*wind**2
    p_drop = (pressure_env - pressure) * 100
    B = rho * e * wind**2 / p_drop
    Vg = (
        np.sqrt(
            ((r / distance) ** B) * (wind**2) * np.exp(1 - (r / distance) ** B)
            + (r**2) * (f**2) / 4
        )
        - (r * f) / 2
    )
    return Vg


start = time.time()

all_winds_path = os.path.join(
    output_dir, "power_intersection", "storm_data", "all_winds", region, sample
)
grid_box = gpd.read_file(
    os.path.join(output_dir, "power_intersection", "regions", f"{region}_unit.gpkg")
)

if len(all_boxes) != 0:
    grid_box = grid_box[
        grid_box["box_id"].isin(all_boxes)
    ]  # filter for boxes that are to be examined only!

totpoints = [len(grid_box)] * len(grid_box)  # for console progress purposes
idx_pts = list(range(len(totpoints)))

sample_num = [sample] * len(grid_box)

print("running wind analysis...")

list_regions = ["NI", "SA", "NA", "EP", "SI", "SP", "WP"]
environmental_pressure = [1006.5, 1014.1, 1014.1, 1008.8, 1010.6, 1008.1, 1008.3]
environ_dict = dict(zip(list_regions, environmental_pressure))

#### load in cyclone tracks for region
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

unique_nh = TC["number_hur"].unique()
TC["sample"] = str(sample)

# TC = TC[TC["number_hur"].isin(nh_func)]  # filter for nhs only

##### change geometry from 0-360 to -180-180
TC["lon"] = TC["lon"].apply(lambda x: x if x <= 180 else x - 360)

hurr_buffer_dist = 1300

distance_arr = np.array([])
unit_path = os.path.join(
    output_dir, "power_intersection", "storm_data", "unit_data", region, sample
)

if not os.path.exists(unit_path):
    os.makedirs(unit_path)


unit_paths_all = (
    dict()
)  # {unit_path: [list of nh in that unit], ... }  Note that this can speed up the loading time for the for loop after next.

for unit in tqdm(grid_box.itertuples(), desc="distances", total=len(grid_box)):

    unique_num = (
        str(unit.longitude)[:7].replace(".", "d").replace("-", "m")
        + "x"
        + str(unit.latitude)[:7].replace(".", "d").replace("-", "m")
        + f"-l{minimum_threshold}c{central_threshold}u{maximum_threshold}"
    )  # implemented such that different units with same ID name (can happen if on second run, there are a different set of all_boxes)
    unit_path_indiv = os.path.join(unit_path, unique_num + ".parquet")
    if not os.path.isfile(unit_path_indiv):
        distance_arr = haversine(unit.longitude, unit.latitude, TC["lon"], TC["lat"])
        TC_sample = TC.copy()
        TC_sample["distance"] = distance_arr
        TC_sample["box_id"] = unit.box_id
        TC_sample["ID_point"] = unit.ID_point
        TC_sample = TC_sample[TC_sample["distance"] <= hurr_buffer_dist]

        nh_sample_unique = list(TC_sample["number_hur"].unique())

        TC_sample.to_parquet(unit_path_indiv, compression="snappy", index=False)

    else:
        print(f"{unit_path_indiv} exists already")

        if nh_split < 25:
            # Option 1: time consuming (in this for loop) but useful if storm_batches is small (will rule out many nh options)
            TC_sample = pd.read_parquet(unit_path_indiv)
            nh_sample_unique = list(TC_sample["number_hur"].unique())

        else:
            # Option 2: quicker (in this for loop) and useful if storm_batches is high (most likely an overlap so less to rule out)
            nh_sample_unique = None

    unit_paths_all.update({unit_path_indiv: nh_sample_unique})


unique_nh_splitlst = [
    unique_nh[i * nh_split : (i + 1) * nh_split]
    for i in range(0, int(len(unique_nh) / nh_split) + 1)
]  # split into lists of length (max) nh_split (is list of lists)
if (
    len(unique_nh_splitlst[-1]) == 0
):  # last one can be [] is nh_split == len(unique_nh_splitlst)
    unique_nh_splitlst = unique_nh_splitlst[:-1]  # remove []


for nh_lst in tqdm(
    unique_nh_splitlst, desc="Storm Damages", total=len(unique_nh_splitlst)
):

    if wind_rerun == False:
        if False in [
            os.path.isfile(
                os.path.join(
                    output_dir,
                    "power_intersection",
                    "storm_data",
                    "all_winds",
                    region,
                    sample,
                    f"TC_r{region}_s{sample}_n{nh}.csv",
                )
            )
            for nh in nh_lst
        ]:
            print("")
        else:
            print("skipping, saved all nh already")
            continue  # skip, all in nh_lst saved already

    TC_all_lst = []
    ss = time.time()
    # for unit_path_indiv in tqdm(unit_paths_all, desc=f'iterating through unit paths for {nh}',total=len(unit_paths_all)):
    for unit_path_indiv, nh_unique in tqdm(
        unit_paths_all.items(), total=len(unit_paths_all), desc="Loading units"
    ):
        if (
            nh_unique != None
        ):  # if it is None, continue because it is unknown which nh are in which units
            if not any([x == y for x in nh_unique for y in nh_lst]):  # if no overlap
                # print("skipping")
                continue  # then dont continue because no point loading as will be empty

        TC_add = pd.read_parquet(unit_path_indiv)
        TC_add = TC_add[TC_add["number_hur"].isin(nh_lst)]
        TC_all_lst.append(TC_add)
    print(f"Time for grid loading: {round((time.time()-ss)/60,3)} mins")

    print("concatenating")
    if len(TC_all_lst) != 0:
        TC_all = pd.concat(TC_all_lst, ignore_index=True)
    else:
        TC_all = []  # dummy

    if len(TC_all) != 0:
        s = time.time()

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
        TC_all = TC_all[
            TC_all["wind_max"] > min_windmax
        ]  ### only with a maximum wind of more than

        ### merge and remove based on mimum distance set
        TC_all = TC_all.merge(min_distance, on="number_hur")
        # TC_all = TC_all[TC_all['distance_min']<max_distance]

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
        TC_all = TC_all[TC_all["wind_location_max"] > min_windlocmax]

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
            ]
        ]
        print("finalising")
        TC_all = TC_all.merge(above20ms, on="number_hur", how="outer").replace(
            np.nan, 0
        )
        TC_all = TC_all.merge(above15ms, on="number_hur", how="outer").replace(
            np.nan, 0
        )
        TC_all["basin"] = region

        print(f"Time for grid processing: {round((time.time()-s)/60,3)} mins")

        if not os.path.exists(all_winds_path):
            os.makedirs(all_winds_path)

        for nh, csv_nh in TC_all.groupby("number_hur"):  #
            print(f"saving {nh}")
            p = os.path.join(all_winds_path, f"TC_r{region}_s{sample}_n{nh}.csv")
            csv_nh.to_csv(p, index=False)

    else:
        print(f"{nh_lst} do not have sufficient unit damage, skipping")

# it is not known in advance how many files will be created by this script
# so, to indicate to snakemake that execution has completed, create a flag file
with open(os.path.join(all_winds_path, "completed.txt"), "w") as fp:
    fp.writelines(f"Winds generated.")

print(f"Total time {round((time.time()-start)/60,3)}")
