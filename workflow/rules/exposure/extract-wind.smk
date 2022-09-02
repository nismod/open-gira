"""
Extracts winds data for specified region sample and year.
"""

import os

try:
    STORM_BATCH_SIZE = int(config["storm_batches"])
    assert STORM_BATCH_SIZE > 0
except:
    raise RuntimeError("storm_batches incorrectly specified in config.yaml file")

WIND_RERUN_BOOL = config["wind_rerun"]
assert WIND_RERUN_BOOL in [True, False]


checkpoint intersect_winds_indiv:
    """Find the .csv files for the wind speed details at each unit.
    IMPORTANT: to reduce computational time, this rule is executed only once and the .py file works out what needs to
               still be calculated. THe output of this rule is limited to rsn_req because when snakemake runs the rule
    it clears all existing files matching the output."""
    input:
        os.path.join(
            config["output_dir"], "power_processed", "world_boxes_metadata.json"
        ),
        os.path.join(
            config["output_dir"], "power_intersection", "regions", "{region}_unit.gpkg"
        ),
        os.path.join(
            config["output_dir"],
            "input",
            "stormtracks",
            "fixed",
            "STORM_FIXED_RETURN_PERIODS_{region}.nc",
        ),
        STORMS_EVENTS,
    params:
        region="{region}",
        sample="{sample}",
        all_boxes_compute=ALL_BOXES,
        memory_storm_split=STORM_BATCH_SIZE,
        wind_rerun=WIND_RERUN_BOOL,
        output_dir=config["output_dir"],
        storm_model_type=config["storm_model_type"],
        wind_file_start=WIND_FILE_START,
        wind_file_end=WIND_FILE_END,
        central_threshold=config["central_threshold"],
        minimum_threshold=config["minimum_threshold"],
        maximum_threshold=config["maximum_threshold"],
    output:
        directory(
            os.path.join(
                config["output_dir"],
                "power_intersection",
                "storm_data",
                "all_winds",
                "{region}",
                "{sample}",
            )
        ),
        os.path.join(
            config["output_dir"],
            "power_intersection",
            "storm_data",
            "all_winds",
            "{region}",
            "{sample}",
            "{region}_{sample}_completed.txt",
        ),
    script:
        os.path.join(
            "..", "..", "scripts", "intersect", "intersect_3_wind_extracter.py"
        )
