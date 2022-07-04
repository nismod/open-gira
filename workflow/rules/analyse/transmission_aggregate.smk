"""Gathers and aggregates statistics on frequency of damage of transmission lines

"""
import os

transmission_out = [
    os.path.join(
        config["output_dir"],
        "power_output",
        "statistics",
        "aggregate",
        "transmission_line_frequency_hit.gpkg",
    ),
    os.path.join(
        config["output_dir"],
        "power_output",
        "statistics",
        "aggregate",
        "transmission_line_reconstruction_costs.gpkg",
    ),
]


rule analyse_transmission:
    input:
        os.path.join(
            stat_path, f"combined_storm_statistics_{config['central_threshold']}.csv"
        ),
    params:
        output_dir=config["output_dir"],
        aggregate_level=config["aggregate_level"],
        region_eval=REGIONS,
        sample_eval=SAMPLES,
        nh_eval=STORMS,
        central_threshold=config["central_threshold"],
    output:
        transmission_out,
    script:
        os.path.join("..", "..", "scripts", "analyse", "transmission_aggregate.py")
