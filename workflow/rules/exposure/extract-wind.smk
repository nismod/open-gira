"""
Extracts winds data for specified region sample and year.
"""

checkpoint intersect_winds_indiv:
    """
    Find windspeeds for every storm for each grid cell.
    """
    conda: "../../../environment.yml"
    input:
        metadata = rules.world_splitter.output.global_metadata,
        unit_grid = "{OUTPUT_DIR}/power_intersection/regions/{STORM_BASIN}_unit.gpkg",
        return_period_map = "{OUTPUT_DIR}/input/storm-ibtracs/fixed/STORM_FIXED_RETURN_PERIODS_{STORM_BASIN}.nc",
        storm_files = STORM_EVENTS,
    params:
        region = lambda wildcards: wildcards.STORM_BASIN,
        sample = lambda wildcards: wildcards.SAMPLE,
        all_boxes_compute = ALL_BOXES,
        memory_storm_split = STORM_BATCH_SIZE,
        wind_rerun = WIND_RERUN_BOOL,
        output_dir = lambda wildcards: wildcards.OUTPUT_DIR,
        storm_file = get_storm_file,
        # wind speed thresholds
        central_threshold = config["central_threshold"],
        minimum_threshold = config["minimum_threshold"],
        maximum_threshold = config["maximum_threshold"],
    output:
        gridded_wind_speeds = directory("{OUTPUT_DIR}/power_intersection/storm_data/all_winds/{STORM_MODEL}/{STORM_BASIN}/{SAMPLE}"),
        completion_flag = "{OUTPUT_DIR}/power_intersection/storm_data/all_winds/{STORM_MODEL}/{STORM_BASIN}/{SAMPLE}/completed.txt",
    script:
        os.path.join(
            "..", "..", "scripts", "intersect", "intersect_3_wind_extracter.py"
        )
