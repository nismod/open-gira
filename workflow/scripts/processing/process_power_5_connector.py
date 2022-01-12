import os
import sys

import geopandas as gpd
import pandas as pd
import json
from tqdm.notebook import tqdm
from importing_modules import *
from process_power_functions import adjbox
import glob
import time


# TODO: remove below lines once testing complete and solely on linux
if "linux" not in sys.platform:
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)

network_paths = glob.glob(os.path.join("data", 'processed', 'all_boxes',"*", "network_*.gpkg"))  # finds the network_{box_id}.gpkg files that exist in all_boxes

for network_path in network_paths:
    s = time.time()
    # extract box numbers
    filename = network_path[network_path.find('network_box_'):]
    box_id = filename[filename.find('box_'):filename.find('.gpkg')]
    id = int(box_id[4:])
    examine = adjbox(id)

    portal_dict = {}
    gdf_idx = gpd.read_file(network_path, layer='nodes')
    if list(gdf_idx.id) != [None]:  # check not dummy

        for id_ex in examine:
            path_test = os.path.join("data", 'processed', 'all_boxes',f"box_{id_ex}", f"network_box_{id_ex}.gpkg")
            if os.path.exists(path_test):
                gdf_ex = gpd.read_file(path_test, layer='nodes')
                if list(gdf_ex.id) != [None]:  # check not dummy
                    # points_ex = gdf_ex.overlay(gdf_idx, how='intersection')  # find overlaps
                    # points_idx = gdf_idx.overlay(gdf_ex, how='intersection')  # find overlaps
                    duplicates = gdf_ex.merge(gdf_idx, left_on='geometry', right_on='geometry', suffixes=('_ex', '_idx'))
                    portal_dict.update({i: j for i, j in zip(duplicates.id_idx, duplicates.id_ex)})  # add {examined_box_point:other_box_point}

    with open(os.path.join("data", 'processed', 'all_boxes', box_id, f"connector_{box_id}.txt"), 'w') as file_ex:
        json.dump(portal_dict, file_ex)

    #print(time.time()-s)

