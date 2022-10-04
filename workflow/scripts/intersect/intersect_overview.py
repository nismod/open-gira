"""Collects all individual (per region, sample) storm json statistics and merges into one csv"""

import pandas as pd
import os

try:
    inputs = snakemake.input  # type: ignore
    output = snakemake.output  # type: ignore
except:  # if the user wishes to see the (partial) statistics before the full analysis, or wishes to view a semi-complete analysis, they can run this file through python3 in the command line. Please note that stat_merger_individual.py MUST be run (once) first
    print("USING DIRECTORY LISTINGS TO FIND STATS (please check manually)")
    print("""!! Assuming output_dir is 'results' !!""")
    output = os.path.join(
        "results", "power_output", "statistics", "combined_storm_statistics__manual.csv"
    )
    stat_path = os.path.join("results", "power_output", "statistics")
    inputs = []
    regions = [
        os.path.basename(path) for path in os.listdir(stat_path) if path[-4:] != ".csv"
    ]  # finds the regions (excludes any csv files in folder)
    for region in regions:
        samples = [
            os.path.basename(path)
            for path in os.listdir(os.path.join(stat_path, region))
        ]
        for sample in samples:
            test_file_inp = os.path.join(
                stat_path,
                region,
                sample,
                f"combined_storm_statistics_{region}_{sample}__manual.csv",
            )
            if os.path.isfile(test_file_inp):
                inputs.append(test_file_inp)

    print(f"Len of inputs is {len(inputs)}")


df = pd.DataFrame()

for input in inputs:
    storm_stats = pd.read_csv(input, keep_default_na=False)
    df = df.append(storm_stats, ignore_index=True)

df.to_csv(str(output), index=False)
print(f"Ouput to {output}")

if len(df) == 0:
    print("Merged all, len=0")
else:
    print("Merged all")
