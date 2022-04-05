"""Takes the large csv and separates the files into smaller (quicker-openable) files"""

import pandas as pd
import os

inputs = snakemake.input
REGIONS = snakemake.params['REGIONS']
SAMPLES = snakemake.params['SAMPLES']


def find_nh(region, sample):
    """Returns a list of all hurricane unique identifiers (nh) in the years years_inp for the given region and sample
    Input:
        years_inp: list of years to check within
        region: str
        sample: str
    Output:
        list of unique nh
    """
    filename = os.path.join(
        "data",
        "stormtracks",
        "events",
        f"STORM_DATA_IBTRACS_{region}_1000_YEARS_{sample}.txt",
    )
    # if (
    #     type(years_inp) != list
    # ):  # if years_inp is not a list, make it a list of length one
    #     years_inp = [years_inp]
    try:
        df = pd.read_csv(filename, header=None)
        df.columns = [
            "year",
            "month",
            "number",
            "step",
            "basin",
            "lat",
            "lon",
            "pressure",
            "wind",
            "radius",
            "cat",
            "landfall",
            "dis_land",
        ]
        # years_inp = [int(year) for year in years_inp]  # ensure integer type
        # df = df[
        #     df["year"].astype(int).isin(years_inp)
        # ]  # filter for selected years (years_inp)
        df["nh_nosample"] = (
            df["year"].astype(int).astype(str)
            + "_"
            + df["number"].astype(int).astype(str)
        )  # new column #_# (first # is year, second # is the storm number of that year)
        nums = list(df["nh_nosample"].unique())  # find all unique #_#
        ret_lst = [
            f"{sample}_{num}" for num in nums
        ]  # add sample for #_#_# (sample is now first #)
    except:
        ret_lst = (
            []
        )  # Need to rerun once the csv file exists (after download).
    return ret_lst

empty_csv = pd.DataFrame({'id':[None]})
nh_completed = set()
for region in REGIONS:
    all_winds_path = os.path.join(
        "data", "intersection", "storm_data", "all_winds", region
    )

    for input in inputs:
        output_files = pd.read_csv(input)
        sample = input.split('_')[-1][1:-4]  # get sample

        sample_path = os.path.join(all_winds_path, sample)
        if not os.path.exists(sample_path):
            os.makedirs(sample_path)


        for nh, csv_nh in output_files.groupby("number_hur"):  #
            print(f"saving {nh}")
            p = os.path.join(sample_path, f"TC_r{region}_s{sample}_n{nh}.csv")
            csv_nh.to_csv(p, index=False)
            nh_completed.add(nh)

        # required for snakemake
        nh_remaining = set(find_nh(region, sample)).difference(nh_completed)
        for nh in nh_remaining:
            p = os.path.join(sample_path, f"TC_r{region}_s{sample}_n{nh}.csv")
            empty_csv.to_csv(p, index=False)
