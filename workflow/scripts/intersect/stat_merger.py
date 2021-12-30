"""Collects all individual storm json statistics and merges into one csv"""

import os
import json
import pandas as pd
import glob
import sys


# TODO: remove below lines once testing complete and solely on linux
if "linux" not in sys.platform:
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)

storm_data_path = os.path.join("data", "intersection", "storm_data")
stormfolders = os.listdir(storm_data_path)

for ii, storm in enumerate(stormfolders):
    storm_file = glob.glob(os.path.join(storm_data_path, storm,"*.txt"))
    assert len(storm_file) == 1
    with open(storm_file[0], 'r') as file:
        storm_stats = json.load(file)
        df_toadd = pd.DataFrame(storm_stats)
    if ii == 0:
        df = df_toadd.copy()
    else:
        df = df.append(df_toadd, ignore_index=True)
df.to_csv(os.path.join("data", "intersection", "combined_storm_statistics.csv"), index=False)
