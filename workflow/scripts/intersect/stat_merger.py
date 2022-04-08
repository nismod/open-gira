"""Collects all individual storm json statistics and merges into one csv"""

import json
import pandas as pd

inputs = snakemake.input
output = snakemake.output

df = pd.DataFrame()

for input in inputs:
    storm_stats = pd.read_csv(input)
    df = df.append(storm_stats, ignore_index=True)

df.to_csv(str(output), index=False)

if len(df) == 0:
    print("Merged all, len=0")
else:
    print("Merged all")
