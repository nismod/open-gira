"""Set of functions which are shared between other scripts"""

import os
import pandas as pd


def find_storm_tots(results_folder, samples, region):
    """Finds the number of storms and number of years for samples and region input"""
    tot_years = 0
    tot_storms = 0

    for sample in samples:
        winds_file = os.path.join(
            results_folder,
            "input",
            "STORM",
            "events",
            f"STORM_DATA_IBTRACS_{region}_1000_YEARS_{sample}.txt",
        )
        winds = pd.read_csv(winds_file, keep_default_na=False)
        winds.columns = [
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
        ]  # https://www.nature.com/articles/s41597-020-0381-2.pdf
        winds["number_hur"] = (
            str(sample)
            + "_"
            + winds["year"].astype(int).astype(str)
            + "_"
            + winds["number"].astype(int).astype(str)
        )

        tot_years += len(winds["year"].unique())
        tot_storms += len(winds["number_hur"].unique())

    return tot_storms, tot_years


def find_storm_files(
    file_type, results_folder, region_eval, sample_eval, nh_eval, thrval
):
    """Will return a list of target.gpkg paths for each storm in the selected region and sample and/or storm identifier, the number of total storms (including those with zero damages) and total years
    Note: file_type options are 'targets' or 'edges' (for transmission lines)
    """
    thrval = str(thrval)
    if file_type == "targets":
        file_name = "targets"
    elif file_type == "edges":
        file_name = "edges_affected"
    else:
        raise RuntimeError("Wrong/no file type entered")

    indiv_storm_path = os.path.join(
        results_folder, "power_intersection", "storm_data", "individual_storms"
    )

    if nh_eval == "None":
        nh_eval = None

    storm_totals = 0
    year_total = 0

    samples_add = []
    if region_eval == None:
        region_eval = [os.path.basename(path) for path in os.listdir(indiv_storm_path)]
        if sample_eval == None:
            for region in region_eval:
                samples_add += [
                    {
                        os.path.basename(path)
                        for path in os.listdir(os.path.join(indiv_storm_path, region))
                    }
                ]

            sample_eval = samples_add[0]
            for samples_test in samples_add[1:]:
                sample_eval = sample_eval.union(
                    samples_test
                )  # only keep overlapping samples. Should not remove samples if correctly entered in config.

            if len(sample_eval) != len(samples_add[0]):
                print(
                    f"Not all samples could be evaluated due to overlap. Checking samples: {sample_eval}"
                )

            for region in region_eval:
                # print(storm_totals, find_storm_tot(results_folder, sample_eval, region))
                storm_add, year_add = find_storm_tots(
                    results_folder, sample_eval, region
                )
                storm_totals += storm_add
                year_total += year_add

    target_paths = []
    for region in region_eval:
        for sample in sample_eval:
            sample = str(sample)
            storms = [
                os.path.basename(path)
                for path in os.listdir(os.path.join(indiv_storm_path, region, sample))
            ]
            if nh_eval == None:
                target_paths_add = [
                    os.path.join(
                        indiv_storm_path,
                        region,
                        sample,
                        file,
                        str(thrval),
                        f"{file_name}__storm_r{region}_s{sample}_n{file[6:]}.gpkg",
                    )
                    for file in storms
                    if os.path.isfile(
                        os.path.join(
                            indiv_storm_path,
                            region,
                            sample,
                            file,
                            str(thrval),
                            f"{file_name}__storm_r{region}_s{sample}_n{file[6:]}.gpkg",
                        )
                    )
                ]  # add only if targets gpkg file exists
            else:
                target_paths_add = [
                    os.path.join(
                        indiv_storm_path,
                        region,
                        sample,
                        file,
                        str(thrval),
                        f"{file_name}__storm_r{region}_s{sample}_n{file[6:]}.gpkg",
                    )
                    for file in storms
                    if os.path.isfile(
                        os.path.join(
                            indiv_storm_path,
                            region,
                            sample,
                            file,
                            str(thrval),
                            f"{file_name}__storm_r{region}_s{sample}_n{file[6:]}.gpkg",
                        )
                    )
                    and file[6:] in nh_eval
                ]  # add only if targets gpkg file exists and in nh_eval

            target_paths += target_paths_add

            storm_add, year_add = find_storm_tots(results_folder, sample_eval, region)
            storm_totals += storm_add
            year_total += year_add

    if len(target_paths) == 0:
        print(f"No {file_name} paths...!")

    print(f"{len(target_paths)} targets and {storm_totals} total storms")
    return target_paths, storm_totals, year_total


def avg(string_input):
    """Returns average key for string_input"""
    return string_input + "_avg"


def sm(string_input):
    """Returns sum key for string_input"""
    return string_input + "_sum"


def ae(string_input):
    """Returns anually expected key for string_input"""
    return string_input + "_anually-expected"


def check_srn(region_eval, sample_eval, nh_eval):
    """Performs some checks and fixes"""

    if region_eval == "None":
        region_eval = None
    if region_eval != None:
        assert all(type(x) == str for x in region_eval)

    if sample_eval == "None":
        sample_eval = None
    if sample_eval != None:
        sample_eval = [str(s) if type(s) != str else s for s in sample_eval]
        assert all(type(x) == str for x in sample_eval)

    if nh_eval == "None":
        nh_eval = None
    if nh_eval != None:
        assert all(type(x) == str for x in nh_eval)

    return region_eval, sample_eval, nh_eval


def traprule(lst, spacing):
    """Trapezium rule"""
    if len(lst) == 0:
        return 0
    if len(lst) == 1:
        return lst[0] * spacing  # dummy
    else:
        return 0.5 * spacing * (lst[0] + lst[-1] + sum(lst[1:-1]))
