"""For each target performs parameter analysis/gathering

"""
import os


rule analyse_targets:
    input:
        os.path.join(
            stat_path, f"combined_storm_statistics_{config['central_threshold']}.csv"
        ),
    params:
        output_dir=config["output_dir"],
        metrics_target=metrics_target,
        top_select=config["top_select"],
        increased_severity_sort=config["increased_severity_sort"],
        region_eval=REGIONS,
        sample_eval=SAMPLES,
        nh_eval=STORMS,
        central_threshold=config["central_threshold"],
    output:
        os.path.join(
            stat_path,
            "aggregate",
            f"targets_geo_top{config['top_select']}{increased_severity_sort_bool}percent.gpkg",
        ),
    script:
        os.path.join(
            "..", "..", "scripts", "analyse", "storm_distribution_empirical_geo.py"
        )
