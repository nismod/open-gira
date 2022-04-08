"""Collects all individual storm json statistics for each region sample and merges into one csv"""

import os
import json
from tqdm import tqdm
import pandas as pd

inputs_stats = snakemake.input
#inputs = snakemake.input
output = str(snakemake.output)
# inputs_stats = inputs[:int(len(inputs)/2)]  # first half are the .txt data files
# if len(inputs_stats) > 1 and len(inputs) > 1:
#     if len(inputs_stats)/2 != len(inputs):
#         print(f"inputs_stats length is {len(inputs_stats)} (/2 is {len(inputs_stats)/2}) and inputs length is {len(inputs)}")
#         basename = os.path.basename(output).split('_')
#         region = basename[-2]
#         sample = basename[-1][:-4]  # remove .csv
#         raise RuntimeError(f"This error is often the result of the incomplete running of the checkpoint rule: intersect_winds_indiv. It is recommended to rerun this rule by deleting the folder data/intersection/storm_data/individual_storms/{region}/{sample}.")

df = pd.DataFrame()

for data_stat in tqdm(inputs_stats, desc='Iterating through stats', total=len(inputs_stats)):
    with open(data_stat, "r") as file:
        storm_stats = json.load(file)
        df_toadd = pd.DataFrame(storm_stats)
        df = df.append(df_toadd, ignore_index=True)

if not os.path.exists(os.path.dirname(output)):
    os.makedirs(output)

df.to_csv(output, index=False)

if len(df) == 0:
    print("Merged, len=0")
else:
    print("Merged")
