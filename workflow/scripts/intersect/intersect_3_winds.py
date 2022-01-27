"""Adapted wind speed file from J Verschuur. Processes stormtracks data and returns the wind speed at each grid location."""

#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import netCDF4 as nc4
import sys
import geopandas as gpd
import os
from shapely.geometry import shape, Polygon, Point
import time
from datetime import date
import fiona
import ast
import json
from shapely.geometry import box
from pathos.multiprocessing import ProcessPool, cpu_count
from tqdm import tqdm


# TODO: remove below lines once testing complete and solely on linux
if "linux" not in sys.platform:
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)
    #list_years = [0,1,2,3,4,5,6,7,8,9]
    sample = 0  # smaller sample size for testing
    region = "SP"
    years_input = ["0"]


else:  # linux
    region = sys.argv[1]
    sample = sys.argv[2]
    years_input = ast.literal_eval(sys.argv[3])




windmaxthreshold = 5  # ignore windspeed less than [m/s]
min_windlocmax = 8  # setting parameter J.V.


def t(num, t):
    print(str(num), time.time()-t)


def haversine(lon1, lat1, lon2, lat2):
    # convert degrees to radians
    lon1 = np.deg2rad(lon1)
    lat1 = np.deg2rad(lat1)
    lon2 = np.deg2rad(lon2)
    lat2 = np.deg2rad(lat2)

    # formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    r_e = 6371
    return c * r_e


def holland_wind_field(r,wind,pressure,pressure_env,distance,lat):
    distance = distance*1000
    r = r*1000
    rho = 1.10
    f = np.abs(1.45842300*10**-4*np.sin(lat))
    e = 2.71828182846
    #p_drop = 2*wind**2
    p_drop = (pressure_env-pressure)*100
    B = rho * e * wind**2 / p_drop
    #B = 1.5
    A = r**B
    Vg = np.sqrt(((r/distance)**B) *(wind**2) * np.exp(1-(r/distance)**B)+(r**2)*(f**2)/4 )-(r*f)/2
    return Vg


def TC_analysis(lat, lon, name, box_id, region, sample, years_func, idx, totpoints):
    print(f"{idx}/{totpoints}")

    list_regions = ['NI','SA','NA','EP','SI','SP','WP']
    environmental_pressure = [1006.5,1014.1,1014.1,1008.8,1010.6,1008.1,1008.3]
    environ_df = pd.DataFrame({'region':list_regions,'pressure':environmental_pressure})

    #### load in cyclone tracks for region
    stormfile = os.path.join("data", "stormtracks", "events", "STORM_DATA_IBTRACS_"+region+"_1000_YEARS_"+sample+".txt")
    TC = pd.read_csv(stormfile, header = None)
    TC.columns = ['year', 'month', 'number', 'step','basin','lat','lon','pressure','wind','radius','cat','landfall','dis_land']  # https://www.nature.com/articles/s41597-020-0381-2.pdf

    TC = TC[TC['year'].astype(int).isin(years_func)]  # filter for year only

    ##### change geometry from 0-360 to -180-180
    mask = (TC['lon'] >180.0)
    TC['lon'][mask] = TC['lon']-360.0

    ### add unique number to cyclone data
    TC['number_hur'] = sample+"_"+TC["year"].astype(int).astype(str) +'_'+ TC["number"].astype(int).astype(str)
    TC['sample'] = str(sample)

    #print(name,country,region,sample)
    # if "linux" not in sys.platform:  # print on windows (when testing)
    #     print(str(idx)+"/"+str(totpoints), "--", sample, "-", lat, lon, name, box_id, region)

    ### background environmental pressure
    pressure_env = environ_df[environ_df['region']==str(region)]

    TC_sample = TC.copy()
    TC_sample['distance'] = haversine(lon, lat, TC_sample['lon'], TC_sample['lat'])
    TC_sample['environ_pressure'] = pressure_env.pressure.iloc[0]

    ### get the maximum wind and minimum distance per event
    max_wind = TC_sample.groupby(['number_hur'])['wind'].max().reset_index().rename(columns = {'wind':'wind_max'})
    min_distance = TC_sample.groupby(['number_hur'])['distance'].min().reset_index().rename(columns = {'distance':'distance_min'})

    ### merge and remove based on lowest threshold for severity hurricane
    TC_sample = TC_sample.merge(max_wind, on = 'number_hur')
    TC_sample  = TC_sample[TC_sample['wind_max']>windmaxthreshold] ### only with a maximum wind of more than  # NOTE: commented

    ### merge and remove based on mimum distance set
    TC_sample = TC_sample.merge(min_distance, on = 'number_hur')
    TC_sample  = TC_sample[TC_sample['distance_min']<2000]  # NOTE: commented

    TC_sample['wind_location'] = holland_wind_field(TC_sample['radius'], TC_sample['wind'], TC_sample['pressure'], TC_sample['environ_pressure'], TC_sample['distance'], TC_sample['lat'])

    max_wind_location = TC_sample.groupby(['number_hur'])['wind_location'].max().reset_index().rename(columns = {'wind_location':'wind_location_max'})

    ### merge and remove based on lowest threshold for severity hurricane at location
    TC_sample = TC_sample.merge(max_wind_location, on = 'number_hur')
    TC_sample  = TC_sample[TC_sample['wind_location_max']>min_windlocmax] # NOTE: commented

    above20ms = TC_sample[TC_sample['wind_location']>20][['wind_location','number_hur']].groupby(['number_hur'])['wind_location'].count().reset_index().rename(columns= {'wind_location':'duration_20ms'})
    above15ms = TC_sample[TC_sample['wind_location']>15][['wind_location','number_hur']].groupby(['number_hur'])['wind_location'].count().reset_index().rename(columns= {'wind_location':'duration_15ms'})

    ### extract the maximum wind speed at point only and associated parameters
    # TODO could add other stats here
    TC_sample_maximum = TC_sample[['number_hur', 'wind_location', 'month', 'year','sample','pressure','distance','wind', 'wind_max']].drop_duplicates(subset = ['number_hur'], keep = 'first')  # removed: .sort_values(by = 'wind_location',ascending = False)
    TC_sample_maximum = TC_sample_maximum.merge(above20ms , on = 'number_hur', how = 'outer').replace(np.nan,0)
    TC_sample_maximum = TC_sample_maximum.merge(above15ms , on = 'number_hur', how = 'outer').replace(np.nan,0)
    TC_sample_maximum['basin'] = str(region)
    TC_sample_maximum['box_id'] = str(box_id)
    TC_sample_maximum['ID_point'] = str(name)

    # below 3 lines are useful for testing but too computationally expensive for full data set
    #TC_sample_maximum['max cat'] = [max(TC[TC['number_hur']==nh_]['cat']) for nh_ in TC['number_hur']]
    #TC_sample_maximum['landfall occurred'] = [max(TC[TC['number_hur']==nh_]['landfall']) for nh_ in TC['number_hur']]
    #assert len(TC_sample_maximum['number_hur'].unique()) == len(TC['number_hur'].unique())  # this is to ensure numbering is 0_0, 0_1, 0_2, 0_3 etc with no jumps (filters for df removed thrice above)

    return TC_sample_maximum



if __name__ == '__main__':  # for windows (due to parallel processing)
    nodesuse = max(1,cpu_count()-2)
    if "linux" not in sys.platform:
        nodesuse = 7

    grid_box = gpd.read_file(os.path.join("data", "intersection", "regions", f"{region}_grid.gpkg"))

    totpoints = [len(grid_box)]*len(grid_box)  # for console progress purposes
    idx_pts = list(range(len(totpoints)))

    sample_num = [sample]*len(grid_box)

    assert type(years_input) == list
    years_toprocess = []
    for year in years_input:
        if not os.path.isfile(os.path.join("data","intersection", "storm_data", "all_winds", f"TC_r{region}_s{sample}_y{year}.csv")):  # {region}_s{sample}_y{year}
            years_toprocess.append(int(year))
    print(f"Years to process: {years_toprocess}")

    all_winds_path = os.path.join("data","intersection", "storm_data", "all_winds")
    if len(years_toprocess) != 0:
        years_num = [years_toprocess]*len(grid_box)



        print("running wind analysis...")
        pool = ProcessPool(nodes=nodesuse)

        s = time.time()
        output = pool.map(TC_analysis, grid_box.latitude, grid_box.longitude, grid_box.ID_point, grid_box.box_id, grid_box.region, sample_num, years_num, idx_pts, totpoints)#, grid_code.idx)
        print(f"Time for grid processing: {round((time.time()-s)/60,3)} mins")


        print("finalising")
        output_files = pd.concat(output)


        if not os.path.exists(all_winds_path):
            os.makedirs(all_winds_path)
        for yr, csv_yr in output_files.groupby('year'):  #
            print(f"saving {int(yr)}")
            years_toprocess.remove(int(yr))
            p = os.path.join(all_winds_path, f'TC_r{region}_s{sample}_y{int(yr)}.csv')
            csv_yr.to_csv(p, index=False)

        csv_empty = pd.DataFrame()
        for year in years_toprocess:  # if any unwritten (because eg storm too weak)
            p = os.path.join(all_winds_path, f'TC_r{region}_s{sample}_y{year}.csv')
            csv_empty.to_csv(p, index=False)

    today = date.today()
    for year in years_input:
        with open(os.path.join(all_winds_path, "log", f"__winds_completed_r{region}_s{sample}_y{year}.txt"), 'w+') as wc:
            json.dump({'date':today.strftime("%d/%m/%Y")}, wc)

    #output_files.to_csv(os.path.join("data","intersection", f'TC_c{code}_r{region}_s{sample}.csv'), index=False)

