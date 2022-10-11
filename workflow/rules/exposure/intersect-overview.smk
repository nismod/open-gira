"""
Merges the gdp loss statistics.

Provides an overview for all storms. This is useful for further analysis
decisions and required to bring together the snakemake wildcards.
"""

import os


def aggregate_input(wildcards):
    # need this to signify that requires checkpoint data
    checkpoint_output = checkpoints.intersect_winds_indiv.get(**wildcards).output.gridded_wind_speeds
    print("checkpoint output:", checkpoint_output)  # in this case: "data"

    # "data2/new_{num}.txt" is the last output before the aggregating file
    # os.path.join(checkpoint_output, "{num}.txt") is the output of the rule that makes all the new wildcards
    # glob_wildcards(os.path.join(checkpoint_output, "{num}.txt")).num  # extracts the actual wildcards (here: {num})
    # expand(...,...) to create the output of the last output before the aggregating file

    return expand(
        os.path.join(
            config["output_dir"],
            "power_intersection",
            "storm_data",
            "individual_storms",
            "{REGION}",
            "{SAMPLE}",
            "storm_{STORM_ID}",
            "{WIND_SPEED_THRESHOLD}",
            "storm_r{REGION}_s{SAMPLE}_n{STORM_ID}.txt",
        ),
        STORM_ID=glob_wildcards(
            os.path.join(checkpoint_output, "TC_r{REGION}_s{SAMPLE}_n{STORM_ID}.csv")
        ).STORM_ID,
        REGION=wildcards.REGION,
        SAMPLE=wildcards.SAMPLE,
        WIND_SPEED_THRESHOLD=wildcards.WIND_SPEED_THRESHOLD,
    )


rule merge_overview_indiv_stats:
    """Use this rule for the combined stats file for each region sample"""
    conda: "../../../environment.yml"
    input:
        aggregate_input,
    output:
        os.path.join(
            "{OUTPUT_DIR}",
            "power_output",
            "statistics",
            "{REGION}",
            "{SAMPLE}",
            "combined_storm_statistics_{REGION}_{SAMPLE}_{WIND_SPEED_THRESHOLD}.csv",
        ),
    script:
        os.path.join(
            "..", "..", "scripts", "intersect", "intersect_overview_individual.py"
        )


rule merge_overview_all_stats:
    """Use this rule for the combined stats file (all)"""
    conda: "../../../environment.yml"
    input:
        STORM_STATS_BY_REGION_SAMPLE_THRESHOLD,
    output:
        STORM_STATS_BY_THRESHOLD,
    params:
        thresholds=WIND_SPEED_THRESHOLDS_MS,
    script:
        os.path.join("..", "..", "scripts", "intersect", "intersect_overview.py")
