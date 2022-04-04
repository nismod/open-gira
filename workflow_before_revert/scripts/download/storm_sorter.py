"""Goes through all storm sets and identifies unique storm itentifiers."""


import sys
import pandas as pd
import os

if 'linux' not in sys.platform:  # TODO

    path = """C:\\Users\\maxor\\Documents\\PYTHON\\GIT\\open-gira"""
    os.chdir(path)


try:
    REGIONS = snakemake.params["REGIONS"]
    SAMPLES = snakemake.params["SAMPLES"]
    YEARS = snakemake.params["YEARS"]
except:
    print("RUNNING MANUAL")
    REGIONS = ["NA", "SP"]
    SAMPLES = [0, 1]
    YEARS = [0, 1, 2, 3]


def find_nh(years_inp, region, sample):
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
    if (
        type(years_inp) != list
    ):  # if years_inp is not a list, make it a list of length one
        years_inp = [years_inp]
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
        years_inp = [int(year) for year in years_inp]  # ensure integer type
        df = df[
            df["year"].astype(int).isin(years_inp)
        ]  # filter for selected years (years_inp)
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
        )  # Need to rerun once the csv file exists (after download). TODO check snakemake checkpoints
    return ret_lst


folder_nh = os.path.join("data", "intersection", "storm_data", "storm_identifiers")
if not os.path.exists(folder_nh):
    os.makedirs(folder_nh)

for region_key in REGIONS:
    for sample_key in SAMPLES:
        nh_to_find = find_nh(YEARS, region_key, sample_key)
        for nh in nh_to_find:
            storm_file = os.path.join(folder_nh, f"r{region_key}_s{sample_key}_n{nh}.txt")
            if not os.path.isfile(storm_file):
                with open(storm_file, 'w') as write_file:
                    write_file.write('')
                print(f"Writing for {region_key} {sample_key}: {storm_file}")
            else:
                print(f"Already written {region_key} {sample_key}")

with open(os.path.join(folder_nh, "dummy.txt"), 'w') as dummy_file:
    dummy_file.write('')