import os

import pandas as pd

def find_storm_tot(results_folder, samples, region):
    for sample in samples:
        winds_file = os.path.join(results_folder, "input", "stormtracks", "events", f"STORM_DATA_IBTRACS_{region}_1000_YEARS_{sample}.txt")
        winds = pd.read_csv(winds_file)
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
        winds['number_hur'] = str(sample) + "_" + winds["year"].astype(int).astype(str) + "_" + winds["number"].astype(int).astype(str)

        return len(winds['number_hur'].unique())


def find_targets(results_folder, region_eval, sample_eval, nh_eval):
    """Will return a list of target.gpkg paths for each storm in the selected region and sample and/or storm identifier and also the number of total storms (including those with zero damages)"""
    indiv_storm_path = os.path.join(results_folder, "power_intersection", "storm_data", "individual_storms")  # TODO change output_dir

    storm_totals = 0

    samples_add = []
    if region_eval == None:
        region_eval = [os.path.basename(path) for path in os.listdir(indiv_storm_path)]
        if sample_eval == None:
            for region in region_eval:
                samples_add += [{os.path.basename(path) for path in os.listdir(os.path.join(indiv_storm_path, region))}]

                  # update




            sample_eval = samples_add[0]
            for samples_test in samples_add[1:]:
                sample_eval = sample_eval.union(samples_test)  # only keep overlapping samples. Should not remove samples if correctly entered in config.

            if len(sample_eval) != len(samples_add[0]):
                print(f'Not all samples could be evaluated due to overlap. Checking samples: {sample_eval}')

            for region in region_eval:
                #print(storm_totals, find_storm_tot(results_folder, sample_eval, region))
                storm_totals += find_storm_tot(results_folder, sample_eval, region)

    target_paths = []
    for region in region_eval:
        for sample in sample_eval:
            storms = [os.path.basename(path) for path in os.listdir(os.path.join(indiv_storm_path, region, sample))]
            if nh_eval == None:
                target_paths_add = [os.path.join(indiv_storm_path, region, sample, file, f"targets__storm_r{region}_s{sample}_n{file[6:]}.gpkg") for file in storms if os.path.isfile(os.path.join(indiv_storm_path, region, sample, file, f"targets__storm_r{region}_s{sample}_n{file[6:]}.gpkg"))]  # add only if targets gpkg file exists
            else:
                target_paths_add = [os.path.join(indiv_storm_path, region, sample, file, f"targets__storm_r{region}_s{sample}_n{file[6:]}.gpkg") for file in storms if os.path.isfile(os.path.join(indiv_storm_path, region, sample, file, f"targets__storm_r{region}_s{sample}_n{file[6:]}.gpkg")) and file[6:] in nh_eval]  # add only if targets gpkg file exists and in nh_eval

            target_paths += target_paths_add

            storm_totals += find_storm_tot(results_folder, [sample], region)  # update


    if len(target_paths) == 0:
        print('No target paths...!')

    print(f'{len(target_paths)} targets and {storm_totals} total storms')
    return target_paths, storm_totals


def avg(string_input):
    """Returns average key for string_input"""
    return string_input+"_avg"


def sm(string_input):
    """Returns sum key for string_input"""
    return string_input+"_sum"