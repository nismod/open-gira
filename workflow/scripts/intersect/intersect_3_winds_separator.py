"""Takes the large csv and separates the files into smaller (quicker-openable) files"""

import pandas as pd
import os

inputs = snakemake.input
REGIONS = snakemake.params['REGIONS']
SAMPLES = snakemake.params['SAMPLES']

print('inputs: ',inputs)

for region in REGIONS:
    all_winds_path = os.path.join(
        "data", "intersection", "storm_data", "all_winds", region
    )

    for input in inputs:
        output_files = pd.read_csv(input)
        sample = input.split('_')[-2][1:]  # get sample

        for nh, csv_nh in output_files.groupby("number_hur"):  #
            print(f"saving {nh}")
            p = os.path.join(all_winds_path, f"TC_r{region}_s{sample}_n{nh}.csv")
            csv_nh.to_csv(p, index=False)