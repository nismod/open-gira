"""
For each target performs parameter analysis/gathering
"""

import os


rule analyse_targets:
    conda: "../../../environment.yml"
    input:
        os.path.join(
            STORM_IMPACT_STATISTICS_DIR, f"combined_storm_statistics_{config['central_threshold']}.csv"
        ),
    params:
        output_dir=config["output_dir"],
        metrics_target=TARGET_ANALYSIS_METRICS,
        top_select=config["top_select"],
        increased_severity_sort=config["increased_severity_sort"],
        region_eval=STORM_BASINS,
        sample_eval=SAMPLES,
        nh_eval=STORMS,
        central_threshold=config["central_threshold"],
    output:
        os.path.join(
            STORM_IMPACT_STATISTICS_DIR,
            "aggregate",
            f"targets_geo_top{config['top_select']}{SORT_BY_INCREASING_SEVERITY}percent.gpkg",
        ),
    script:
        os.path.join(
            "..", "..", "scripts", "analyse", "storm_distribution_empirical_geo.py"
        )
