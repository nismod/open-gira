"""Collects all individual storm json statistics for each region sample and merges into one csv"""

import os
import json
import pandas as pd

inputs = snakemake.inputs
output = str(snakemake.outputs)
inputs_stats = inputs[:int(len(inputs)/2)]  # first half are the .txt data files
assert len(inputs_stats)/2 == len(inputs)

df = pd.DataFrame()

for data_stat in inputs_stats:
    with open(data_stat, "r") as file:
        storm_stats = json.load(file)
        df_toadd = pd.DataFrame(storm_stats)
        df = df.append(df_toadd, ignore_index=True)

if not os.path.exists(os.path.dirname(output)):
    os.makedirs(output)

df.to_csv(output)

if len(df) == 0:
    print("Merged, len=0")
else:
    print("Merged")
