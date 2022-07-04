"""Based on percentile from config file, will copy to percentile folder in statistics

"""
import os

percentile_out = directory(
    os.path.join(config["output_dir"], "power_output", "statistics", "percentile")
)


rule analyse_percentile:
    input:
        stat_csv,
    params:
        output_dir=config["output_dir"],
        region_eval=REGIONS,
        sample_eval=SAMPLES,
        nh_eval=STORMS,
        metrics_target=metrics_target,
        central_threshold=config["central_threshold"],
        percentile=config["percentile"],
    output:
        percentile_out,
    script:
        os.path.join("..", "..", "scripts", "analyse", "select_percentile.py")
