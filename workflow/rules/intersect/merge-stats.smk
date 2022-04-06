"""Merges the gdp loss statistics.

"""
import os

stat_csv = os.path.join("data", "intersection", "statistics", "combined_storm_statistics.csv")


def aggregate_input(wildcards):
    checkpoint_output = checkpoints.intersect_winds_indiv.get(**wildcards).output[0] #need this to signify that requires checkpoint data
    print('checkpoint output:', checkpoint_output)  # in this case: "data"

    # "data2/new_{num}.txt" is the last output before the aggregating file
    # os.path.join(checkpoint_output, "{num}.txt") is the output of the rule that makes all the new wildcards
    # glob_wildcards(os.path.join(checkpoint_output, "{num}.txt")).num  # extracts the actual wildcards (here: {num})
    # expand(...,...) to create the output of the last output before the aggregating file

    ret = expand([os.path.join(
                "data",
                "intersection",
                "storm_data",
                "individual_storms",
                "{region}",
                "storm_{nh}",
                "storm_r{region}_s{sample}_n{nh}.txt",
            ),
            os.path.join(
                "data",
                "intersection",
                "storm_data",
                "individual_storms",
                "{region}",
                "storm_{nh}",
                "storm_track_r{region}_s{sample}_n{nh}.gpkg",
            )],
        nh=glob_wildcards(os.path.join(checkpoint_output, "TC_r{region}_s{sample}_n{nh}.csv")).nh,
        region=wildcards.region,
        sample=wildcards.sample
    )
    #print(ret)
    return ret
#os.path.join("data", "intersection", "storm_data", "all_winds", "{region}", "{sample}", )



# storm_details_all = expand(
#     [
#         os.path.join(
#             "data",
#             "intersection",
#             "storm_data",
#             "individual_storms",
#             "{region}",
#             "storm_{nh}",
#             "storm_r{region}_s{sample}_n{nh}.txt",
#         ),
#         os.path.join(
#             "data",
#             "intersection",
#             "storm_data",
#             "individual_storms",
#             "{region}",
#             "storm_{nh}",
#             "storm_track_r{region}_s{sample}_n{nh}.gpkg",
#         ),
#     ],
#     region=REGIONS,
#     sample=SAMPLES,
#     nh=find_nh_mult(YEARS, REGIONS, SAMPLES)[2],
# )
#print(storm_details_all)

rule merge_indiv_stats:
    """Use this rule for the combined stats file for each region sample"""
    input:
        aggregate_input,
    output:
        os.path.join("data", "intersection", "statistics", "{region}", "{sample}", "combined_storm_statistics_{region}_{sample}.csv")
    script:
            os.path.join("..", "..", "scripts", "intersect", "stat_merger_individual.py"
        )


rule merge_all_stats:
    """Use this rule for the combined stats file (all)"""
    input:
        expand(os.path.join("data", "intersection", "statistics", "{region}", "{sample}", "combined_storm_statistics_{region}_{sample}.csv"), region=REGIONS, sample=SAMPLES)
    output:
        stat_csv
    script:
            os.path.join("..", "..", "scripts", "intersect", "stat_merger.py"
        )