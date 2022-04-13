"""Extracts winds data for specified region sample and year.

"""
import os

try:
    storm_batch_value = float(config["storm_batches"])
    assert 1 <= storm_batch_value
except:
    raise RuntimeError("storm_batches incorrectly specified in config.yaml file")

wind_rerun_bool = config["wind_rerun"]
assert wind_rerun_bool in [True, False]


checkpoint intersect_winds_indiv:
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
        region="{region}",
        sample="{sample}",
        all_boxes_compute=all_boxes,
        memory_storm_split=storm_batch_value,
        wind_rerun=wind_rerun_bool,
    output:
        directory(
            os.path.join(
                "data",
                "intersection",
                "storm_data",
                "all_winds",
                "{region}",
                "{sample}",
            )
        ),
        os.path.join(
            "data",
            "intersection",
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
