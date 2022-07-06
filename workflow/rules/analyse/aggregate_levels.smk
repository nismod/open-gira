"""Takes a gpkg file and aggregates to chosen level

"""
import os

aggregate_levels_out = os.path.join(
    stat_path,
    "aggregate",
    f"targets_geo_top{config['top_select']}{increased_severity_sort_bool}percent_aggregated_region.gpkg",
)


rule analyse_aggregate_levels:
    input:
        os.path.join(
            stat_path,
            "aggregate",
            f"targets_geo_top{config['top_select']}{increased_severity_sort_bool}percent.gpkg",
        ),
        os.path.join(
            config["output_dir"], "input", "adminboundaries", f"gadm36_levels.gpkg"
        ),
    params:
        output_dir=config["output_dir"],
        metrics_target=metrics_target,
        top_select=config["top_select"],
        increased_severity_sort=config["increased_severity_sort"],
        aggregate_level=config["aggregate_level"],
    output:
        aggregate_levels_out,
    script:
        os.path.join("..", "..", "scripts", "analyse", "storm_aggregate_levels.py")
