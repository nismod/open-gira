"""Collects all individual storm json statistics and merges into one csv"""

import json
import pandas as pd

inputs = snakemake.inputs
output = snakemake.outputs

df = pd.DataFrame()

for input in inputs:
    with open(input, "r") as file:
        storm_stats = json.load(file)
        df_toadd = pd.DataFrame(storm_stats)
        df = df.append(df_toadd, ignore_index=True)

df.to_csv(str(output))

if len(df) == 0:
    print("Merged, len=0")
else:
    print("Merged")
