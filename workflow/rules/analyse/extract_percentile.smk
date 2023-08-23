"""
Based on percentile from config file, will copy to percentile folder in statistics
"""

import os


PERCENTILE_OUT = directory(
    os.path.join(config["output_dir"], "power_output", "statistics", "percentile")
)


rule analyse_percentile:
    input:
        STORM_STATS_BY_THRESHOLD,
    params:
        output_dir=config["output_dir"],
        region_eval=STORM_BASINS,
        sample_eval=SAMPLES,
        nh_eval=STORMS,
        metrics_target=TARGET_ANALYSIS_METRICS,
        central_threshold=config["central_threshold"],
        percentile=config["percentile"],
    output:
        PERCENTILE_OUT,
    script:
        os.path.join("..", "..", "scripts", "analyse", "select_percentile.py")
