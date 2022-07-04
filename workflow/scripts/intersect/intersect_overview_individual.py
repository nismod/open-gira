"""Collects all individual storm json statistics for each region sample and merges into one csv"""





import os
import json
from tqdm import tqdm
import pandas as pd


try:
    inputs_stats = snakemake.input
    outputs = [str(snakemake.output)]  # made list to counter for except: case
except:
    raise NotImplementedError("use snakemake")

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
