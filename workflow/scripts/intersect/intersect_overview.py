"""Collects all individual (per region, sample) storm json statistics and merges into one csv"""

import pandas as pd
import os

try:
    inputs = snakemake.input
    outputs = snakemake.output
    thresholds = snakemake.params["thresholds"]
except:
    raise NotImplementedError("use snakemake")


for thrval in thresholds:
    inputs_sorted = [i for i in inputs if str(thrval) in os.path.basename(i)]
    output_list = [i for i in outputs if str(thrval) in os.path.basename(i)]
    assert len(output_list) == 1
    output = output_list[0]
    df = pd.DataFrame()
    for input in inputs_sorted:
        storm_stats = pd.read_csv(input, keep_default_na=False)
        df = df.append(storm_stats, ignore_index=True)

    df.to_csv(str(output), index=False)
    print(f"Output to {output}")

    if len(df) == 0:
        print("Merged, len=0")
    else:
        print("Merged")

print("Merged all")
