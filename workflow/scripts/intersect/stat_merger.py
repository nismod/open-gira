"""Collects all individual storm json statistics and merges into one csv"""

import os
import json
import pandas as pd
import glob
import sys


if "linux" not in sys.platform:
    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)

storm_data_path = os.path.join("data", "intersection", "storm_data", "individual_storms")
stormfolders = os.listdir(storm_data_path)
stormfolders = [foldername for foldername in stormfolders if foldername[:6]=='storm_']  # keep only storm ones

df = pd.DataFrame()
for ii, storm in enumerate(stormfolders):
    storm_file = glob.glob(os.path.join(storm_data_path, storm,"*.txt"))
    if len(storm_file) == 1:
        with open(storm_file[0], 'r') as file:
            storm_stats = json.load(file)
            df_toadd = pd.DataFrame(storm_stats)
            df = df.append(df_toadd, ignore_index=True)
df.to_csv(os.path.join("data", "intersection", "combined_storm_statistics.csv"), index=False)
if len(df) == 0:
    print('Merged, len=0')
else:
    print('Merged')
