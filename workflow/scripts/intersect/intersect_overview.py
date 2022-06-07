"""Collects all individual (per region, sample) storm json statistics and merges into one csv"""

import pandas as pd
import os

try:
    inputs = snakemake.input
    outputs = snakemake.output
    thresholds = snakemake.params['thresholds']
except:  # if the user wishes to see the (partial) statistics before the full analysis, or wishes to view a semi-complete analysis, they can run this file through python3 in the command line. Please note that stat_merger_individual.py MUST be run (once) first
    raise NotImplementedError("use snakemake")
    # print("USING DIRECTORY LISTINGS TO FIND STATS (please check manually)")
    # thrvalcen = 40
    # thrvalmin = 35
    # thrvalmax = 45
    # thrval_lst = [thrvalcen, thrvalmin, thrvalmax]
    # print(f"""!! Assuming output_dir is 'results' and storm thrvals are {thrvalcen}, {thrvalmin}, {thrvalmax}!!""")
    # output = os.path.join(
    #     "results", "power_output", "statistics", "combined_storm_statistics__manual.csv"
    # )
    # stat_path = os.path.join("results", "power_output", "statistics")
    # inputs = []
    # regions = [
    #     os.path.basename(path) for path in os.listdir(stat_path) if os.path.basename(path) in ["EP", "NA", "NI", "SI", "SP", "WP"]
    # ]  # finds the regions (excludes any csv files in folder)
    # for region in regions:
    #     samples = [
    #         os.path.basename(path)
    #         for path in os.listdir(os.path.join(stat_path, region))
    #     ]
    #     for sample in samples:
    #         for thrval in thrval_lst:
    #             test_file_inp = os.path.join(
    #                 stat_path,
    #                 region,
    #                 sample,
    #                 f"combined_storm_statistics_{region}_{sample}_{thrval}__manual.csv",
    #             )
    #             if os.path.isfile(test_file_inp):
    #                 inputs.append(test_file_inp)
    #
    # print(f"Len of inputs is {len(inputs)}")
    #



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

print('Merged all')
