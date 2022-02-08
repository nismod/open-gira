"""Extracts winds data for specified region sample and year.

"""
import os
import pandas as pd


def required_nh_remaining(rsn):
    """Checks which storms (unique identifier number_hur, abbreviated: nh) needs to be evaluated. Note that each of
       the lists' element locations correspond to the same event.
    Input:
        rsn: list [list of regions, list of samples, list of nh]
    Output:
        remaining elements: list [list of remaining regions, list of remaining samples, list of remaining nh]
    """
    region_all = rsn[0]
    sample_all = rsn[1]
    nh_all = rsn[2]
    nh_completed_files = []
    for region in region_all:
        nh_completed_files += glob(
            os.path.join(
                DATA_DIR, "intersection", "storm_data", "all_winds", region, "*csv"
            )
        )  # find which wind speeds have been completed already
    nh_completed = [
        file[file.find("_n") + 2 : -4] for file in nh_completed_files
    ]  # single out nh identifiers from directory files
    indices_remaining = [
        x for x, nh in enumerate(nh_all) if nh not in nh_completed
    ]  # list all indices that are not already completed
    region_remaining = [
        r for x, r in enumerate(region_all) if x in indices_remaining
    ]  # select based on indices
    sample_remaining = [
        s for x, s in enumerate(sample_all) if x in indices_remaining
    ]  # select based on indices
    nh_remaining = [
        nh for x, nh in enumerate(nh_all) if x in indices_remaining
    ]  # select based on indices
    print("Remaining storms to evaluate wind speeds:", nh_remaining)
    return [region_remaining, sample_remaining, nh_remaining]


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


def find_nh_mult(years_inp, regions_inp, samples_inp):
    """Given multiple regions and multiple samples, find the nh identifiers over the years years_inp
    Input:
        years_inp: list of years
        regions_inp: list of regions
        samples_inp: list of samples
    Output:
        list [list of regions, list of samples, list of nh]
    NOTE: The indices of each list correspond to one another. E.g. the second element of list of nh is in the region
          given by the second element of list of regions and in the sample of the second element of list of samples.
    """
    nh_indiv = []
    region_indiv = []
    sample_indiv = []
    for region in regions_inp:
        for sample in samples_inp:
            add_nh = find_nh(
                years_inp, region, sample
            )  # for each region and sample, find the unique nh identifiers
            nh_indiv += add_nh  # store this
            region_indiv += [region] * len(
                add_nh
            )  # corresponds to the selected region (hence multiply by length to ensure all list elements correspond)
            sample_indiv += [sample] * len(add_nh)  # see above
    return [region_indiv, sample_indiv, nh_indiv]


# TC_all is all outputs expected from the wind speed analysis
TC_all = expand(
    os.path.join(
        "data",
        "intersection",
        "storm_data",
        "all_winds",
        "{region}",
        "TC_r{region}_s{sample}_n{nh}.csv",
    ),
    region=REGIONS,
    sample=SAMPLES,
    nh=find_nh_mult(YEARS, REGIONS, SAMPLES)[2],
)

rsn_req = required_nh_remaining(
    find_nh_mult(YEARS, REGIONS, SAMPLES)
)  # [region list, sample list, nh list] for all nh storms that have NOT had their wind speed calculations for the .csv


rule intersect_winds_indiv:
    """Find the .csv files for the wind speed details at each unit. 
    IMPORTANT: to reduce computational time, this rule is executed only once and the .py file works out what needs to
               still be calculated. THe output of this rule is limited to rsn_req because when snakemake runs the rule
    it clears all existing files matching the output."""
    input:
        os.path.join("data", "processed", "world_boxes_metadata.txt"),
        os.path.join("data", "intersection", "regions", "{region}_unit.gpkg"),
        os.path.join(
            "data", "stormtracks", "fixed", "STORM_FIXED_RETURN_PERIODS_{region}.nc"
        ),
        os.path.join(
            "data",
            "stormtracks",
            "events",
            "STORM_DATA_IBTRACS_{region}_1000_YEARS_{sample}.txt",
        ),
    params:
        nh_compute=lambda wildcards: str(
            find_nh(YEARS, wildcards.region, wildcards.sample)
        ),
    output:
        [
            os.path.join(
                "data",
                "intersection",
                "storm_data",
                "all_winds",
                "{region}" "TC_r{region}_s{sample}_n" + f"{nh}.csv",
            )
            for nh in rsn_req[2]
        ],
    shell:
        (
            "python3 "
            + os.path.join("workflow", "scripts", "intersect", "intersect_3_winds.py")
            + " {wildcards.region} {wildcards.sample} "
            + '"""'
            + "{params.nh_compute}"
            + '"""'
        )


rule intersect_wind:
    """For all elements"""
    input:
        TC_all,
