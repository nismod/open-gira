"""
Takes a gpkg file and aggregates to chosen level
"""

import os

AGGREGATE_LEVELS_OUT = os.path.join(
    STORM_IMPACT_STATISTICS_DIR,
    "aggregate",
    f"targets_geo_top{config['top_select']}{SORT_BY_INCREASING_SEVERITY}percent_aggregated_region.gpkg",
)


rule analyse_aggregate_levels:
    conda: "../../../environment.yml"
    input:
        os.path.join(
            STORM_IMPACT_STATISTICS_DIR,
            "aggregate",
            f"targets_geo_top{config['top_select']}{SORT_BY_INCREASING_SEVERITY}percent.gpkg",
        ),
        os.path.join(
            config["output_dir"], "input", "admin-boundaries", f"gadm36_levels.gpkg"
        ),
    params:
        output_dir=config["output_dir"],
        metrics_target=TARGET_ANALYSIS_METRICS,
        top_select=config["top_select"],
        increased_severity_sort=config["increased_severity_sort"],
        aggregate_level=config["aggregate_level"],
    output:
        AGGREGATE_LEVELS_OUT,
    script:
        os.path.join("..", "..", "scripts", "analyse", "storm_aggregate_levels.py")
