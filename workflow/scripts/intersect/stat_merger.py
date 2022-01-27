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

years_data_path = os.path.join("data", "intersection", "storm_data", "damages")

year_files = glob.glob(os.path.join(years_data_path, "*.txt"))
df = pd.DataFrame()
for year in year_files:
    with open(year, 'r') as file:
        storm_stats = json.load(file)
        for stats_add in list(storm_stats.values()):
            df_toadd = pd.DataFrame(stats_add)
            df = df.append(df_toadd, ignore_index=True)
if len(df) != 0:
    df.to_csv(os.path.join("data", "intersection", "combined_storm_statistics.csv"), index=False)
