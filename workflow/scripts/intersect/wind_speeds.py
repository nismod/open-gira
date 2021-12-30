"""Adapted wind speed file from J Verschuur. Processes stormtracks data and returns the wind speed at each grid location."""

#!/usr/bin/env python
# coding: utf-8


import numpy as np
from shapely.geometry import Point, LineString
from shapely.geometry.polygon import Polygon
from geopandas import GeoSeries
import pandas as pd
import json
import sys
import geopandas as gpd
import os
from shapely.geometry import shape, Polygon
import fiona
from pathos.multiprocessing import ProcessPool, cpu_count

spres = 0.25  # spatial resolution  TODO use return period!

# TODO: remove below lines once testing complete and solely on linux
if "linux" not in sys.platform:
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)
    #list_years = [0,1,2,3,4,5,6,7,8,9]
    sample = 0  # smaller sample size for testing
    code = "PHL"
    region = "WP"
    #storm_num = "0"
else:  # linux
    code = sys.argv[1]
    region = sys.argv[2]
    sample = sys.argv[3]
    #storm_num = sys.argv[4]


def return_lat_lon_bound(geometry_polygon):
    array_coordinates = np.array(geometry_polygon.exterior.coords)
    lon = [array_coordinates[i][0] for i in range(0,len(array_coordinates))]
    lat = [array_coordinates[i][1] for i in range(0,len(array_coordinates))]
    return np.min(lon),np.max(lon),np.min(lat),np.max(lat)


def make_grid_points(country_geodata,res = 0.25):
    country_df = country_geodata.explode(index_parts=True)  # added index_parts=True due to FutureWarning
    lon_min_list = []
    lon_max_list = []
    lat_min_list = []
    lat_max_list = []
    for i in range(0,len(country_df)):
        lon_min,lon_max,lat_min,lat_max = return_lat_lon_bound(country_df.iloc[i].geometry)
        lon_min_list.append(lon_min)
        lon_max_list.append(lon_max)
        lat_min_list.append(lat_min)
        lat_max_list.append(lat_max)
    country_df['lon_min'] =  lon_min_list
    country_df['lon_max'] =  lon_max_list
    country_df['lat_min'] =  lat_min_list
    country_df['lat_max'] =  lat_max_list
    ### find the max and min
    min_lat_c = country_df['lat_min'].min()
    max_lat_c = country_df['lat_max'].max()
    min_lon_c = country_df['lon_min'].min()
    max_lon_c = country_df['lon_max'].max()

    #### get the number of grids, assume a grid cell is around 30km, or around 0.25 degree
    vertical_grid = np.ceil((max_lat_c - min_lat_c)/res)  ### round upwards
    horizontal_grid = np.ceil((max_lon_c - min_lon_c)/res) ### round upwards

    delta_vertical = (max_lat_c - min_lat_c)/vertical_grid
    delta_horizontal = (max_lon_c - min_lon_c)/horizontal_grid

    vertical_centroid = np.linspace(min_lat_c+delta_vertical/2,max_lat_c-delta_vertical/2,vertical_grid.astype(int))
    horizontal_centroid = np.linspace(min_lon_c+delta_horizontal/2,max_lon_c-delta_horizontal/2,horizontal_grid.astype(int))


    point_df = pd.DataFrame()
    ### loop through them and make centroids:
    for i in vertical_centroid:
        for j in horizontal_centroid:
            point_df_add = pd.DataFrame({'longitude':[j],'latitude':[i]})
            point_df = pd.concat([point_df,point_df_add])
    gdf = gpd.GeoDataFrame(point_df, geometry=gpd.points_from_xy(point_df.longitude, point_df.latitude))

    return gdf


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


def TC_analysis(lat, lon, name, country, region, sample):
    list_regions = ['NI','SA','NA','EP','SI','SP','WP']
    environmental_pressure = [1006.5,1014.1,1014.1,1008.8,1010.6,1008.1,1008.3]
    environ_df = pd.DataFrame({'region':list_regions,'pressure':environmental_pressure})

    threshold1 = 10.0
    #### load in cyclone tracks for region
    stormfile = os.path.join("data", "stormtracks", "events", "STORM_DATA_IBTRACS_"+region+"_1000_YEARS_"+sample+".txt")
    TC = pd.read_csv(stormfile, header = None)
    TC.columns = ['year', 'month', 'number', 'step','basin','lat','lon','pressure','wind','radius','cat','landfall','dis_land']

    ##### change geometry from 0-360 to -180-180
    mask = (TC['lon'] >180.0)
    TC['lon'][mask] = TC['lon']-360.0

    ### add unique number to cyclone data
    TC['number_hur'] = TC["year"].astype(int).astype(str) +'_'+ TC["number"].astype(int).astype(str)
    TC['sample'] = str(sample)

    #print(name,country,region,sample)
    if "linux" not in sys.platform:  # print on windows (when testing)
        print(sample, "--", lat, lon, name, country, region)

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
    #TC_sample  = TC_sample[TC_sample['wind_max']>threshold1] ### only with a maximum wind of more than  # NOTE: commented

    ### merge and remove based on mimum distance set
    TC_sample = TC_sample.merge(min_distance, on = 'number_hur')
    #TC_sample  = TC_sample[TC_sample['distance_min']<2000]  # NOTE: commented

    TC_sample['wind_location'] = holland_wind_field(TC_sample['radius'], TC_sample['wind'], TC_sample['pressure'], TC_sample['environ_pressure'], TC_sample['distance'], TC_sample['lat'])

    max_wind_location = TC_sample.groupby(['number_hur'])['wind_location'].max().reset_index().rename(columns = {'wind_location':'wind_location_max'})

    ### merge and remove based on lowest threshold for severity hurricane at location
    TC_sample = TC_sample.merge(max_wind_location, on = 'number_hur')
    #TC_sample  = TC_sample[TC_sample['wind_location_max']>8] # NOTE: commented

    above20ms = TC_sample[TC_sample['wind_location']>20][['wind_location','number_hur']].groupby(['number_hur'])['wind_location'].count().reset_index().rename(columns= {'wind_location':'duration_20ms'})
    above15ms = TC_sample[TC_sample['wind_location']>15][['wind_location','number_hur']].groupby(['number_hur'])['wind_location'].count().reset_index().rename(columns= {'wind_location':'duration_15ms'})

    ### extract the maximum wind speed at point only and associated parameters
    TC_sample_maximum = TC_sample[['number_hur', 'wind_location', 'month', 'year','sample','pressure','distance','wind', 'wind_max']].sort_values(by = 'wind_location',ascending = False).drop_duplicates(subset = ['number_hur'], keep = 'first')
    TC_sample_maximum = TC_sample_maximum.merge(above20ms , on = 'number_hur', how = 'outer').replace(np.nan,0)
    TC_sample_maximum = TC_sample_maximum.merge(above15ms , on = 'number_hur', how = 'outer').replace(np.nan,0)
    TC_sample_maximum['basin'] = str(region)
    TC_sample_maximum['country'] = str(country)
    TC_sample_maximum['ID_point'] = str(name)

    assert len(TC_sample_maximum['number_hur'].unique()) == len(TC['number_hur'].unique())  # this is to ensure numbering is 0_0, 0_1, 0_2, 0_3 etc with no jumps (filters for df removed thrice above)

    return TC_sample_maximum



if __name__ == '__main__':  # for windows (due to parallel processing)
    print("opening files")
    with fiona.open(os.path.join("data","adminboundaries", f"gadm36_{code}.gpkg"), "r", layer=3) as src:
        code_geoms = []
        for feature in src:
            code_geoms.append(shape(feature['geometry']))
        country = gpd.GeoDataFrame({'geometry':code_geoms})

    ### create an artifical grid
    grid_code = make_grid_points(country, spres)

    grid_code = grid_code.overlay(country, how='intersection')

    grid_code = grid_code.reset_index(drop = True)
    grid_code['ID_point'] = grid_code.index

    ### extract points inside
    grid_code['country'] = code
    grid_code.to_file(os.path.join("data","intersection", f"grid_{code}.gpkg"),driver = 'GPKG')
    grid_code['region'] = region

    grid_code['sample_num'] = str(sample)

    #grid_code['idx'] = storm_num

    nodesuse = max(1,cpu_count()-2)
    if "linux" not in sys.platform:
        nodesuse = 3

    print("running wind analysis...")
    pool = ProcessPool(nodes=nodesuse)
    output = pool.map(TC_analysis, grid_code.latitude, grid_code.longitude, grid_code.ID_point, grid_code.country, grid_code.region, grid_code.sample_num)#, grid_code.idx)
    print("finalising")
    output_files = pd.concat(output)
    output_files.to_csv(os.path.join("data","intersection", f'TC_c{code}_r{region}_s{sample}.csv'), index=False)

