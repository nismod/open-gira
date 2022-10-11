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
    """
    Find windspeeds for every storm for each grid cell.
    """
    conda: "../../../environment.yml"
    input:
        metadata = rules.world_splitter.output.global_metadata,
        unit_grid = "{OUTPUT_DIR}/power_intersection/regions/{REGION}_unit.gpkg",
        return_period_map = "{OUTPUT_DIR}/input/stormtracks/fixed/STORM_FIXED_RETURN_PERIODS_{REGION}.nc",
        storm_files = STORMS_EVENTS,
    params:
        region = lambda wildcards: wildcards.REGION,
        sample = lambda wildcards: wildcards.SAMPLE,
        all_boxes_compute = ALL_BOXES,
        memory_storm_split = STORM_BATCH_SIZE,
        wind_rerun = WIND_RERUN_BOOL,
        output_dir = lambda wildcards: wildcards.OUTPUT_DIR,
        storm_model_type = STORM_MODEL,
        wind_file_start = WIND_FILE_START,
        wind_file_end = WIND_FILE_END,
        # wind speed thresholds
        central_threshold = config["central_threshold"],
        minimum_threshold = config["minimum_threshold"],
        maximum_threshold = config["maximum_threshold"],
    output:
        gridded_wind_speeds = directory("{OUTPUT_DIR}/power_intersection/storm_data/all_winds/{REGION}/{SAMPLE}"),
        completion_flag = "{OUTPUT_DIR}/power_intersection/storm_data/all_winds/{REGION}/{SAMPLE}/completed.txt",
    script:
        os.path.join(
            "..", "..", "scripts", "intersect", "intersect_3_wind_extracter.py"
        )
