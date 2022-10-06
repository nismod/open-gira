"""Collects all individual storm json statistics for each region sample and merges into one
csv
"""
import json
import os

import pandas as pd
from tqdm import tqdm

try:
    inputs_stats = snakemake.input  # type: ignore
    outputs = [str(snakemake.output)]  # type: ignore
except:
    raise NotImplementedError("use snakemake")

# outputs made list to counter for except: case
for (
    output
) in (
    outputs
):  # loop to counter for non-snakemake runs (in which all need to be executed during script as no wildcards possible)
    df = pd.DataFrame()

    for data_stat in tqdm(
        inputs_stats, desc="Iterating through stats", total=len(inputs_stats)
    ):
        with open(data_stat, "r") as file:
            storm_stats = json.load(file)
            df_toadd = pd.DataFrame(storm_stats)

            df = df.append(df_toadd, ignore_index=True)

    if not os.path.exists(os.path.dirname(output)):
        os.makedirs(os.path.dirname(output))

    df.to_csv(output, index=False)
    print(f"Output to {output}")

    if len(df) == 0:
        print("Merged, len=0")
    else:
        print("Merged")
