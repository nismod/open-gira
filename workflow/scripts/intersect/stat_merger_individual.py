"""Collects all individual storm json statistics for each region sample and merges into one csv"""

import os
import json
from tqdm import tqdm
import pandas as pd

try:
    inputs_stats = snakemake.input
    outputs = [str(snakemake.output)]  # made list to counter for except: case
except:  # if the user wishes to see the (partial) statistics before the full analysis, or wishes to view a semi-complete analysis, they can run this file through python3 in the command line. Please note that stat_merger.py must be run (once) after this to see ALL collated storm results.
    print("USING DIRECTORY LISTINGS TO FIND STATS (please check manually)")
    print("""!! Assuming output_dir is 'results' !!""")
    indiv_storm_path = os.path.join(
        "results", "power_intersection", "storm_data", "individual_storms"
    )
    inputs_stats = []
    regions = [os.path.basename(path) for path in os.listdir(indiv_storm_path)]
    reg_sam = []
    for region in regions:
        samples = [
            os.path.basename(path)
            for path in os.listdir(os.path.join(indiv_storm_path, region))
        ]
        for sample in samples:
            storms = [
                os.path.basename(path)
                for path in os.listdir(os.path.join(indiv_storm_path, region, sample))
            ]
            new_stats = [
                os.path.join(
                    indiv_storm_path,
                    region,
                    sample,
                    file,
                    f"storm_r{region}_s{sample}_n{file[6:]}.txt",
                )
                for file in storms
                if os.path.isfile(
                    os.path.join(
                        indiv_storm_path,
                        region,
                        sample,
                        file,
                        f"storm_r{region}_s{sample}_n{file[6:]}.txt",
                    )
                )
            ]  # add only if txt file exists
            if len(new_stats) >= 1:
                reg_sam += [(region, sample)]
            inputs_stats += new_stats
    print(f"Len of inputs_stats is {len(inputs_stats)}")
    outputs = [
        os.path.join(
            "results",
            "power_output"
            "statistics",
            f"{region}",
            f"{sample}",
            f"combined_storm_statistics_{region}_{sample}__manual.csv",
        )
        for region, sample in reg_sam
    ]  # list of outputs
    print(f"Len of outputs is {len(outputs)}")


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
