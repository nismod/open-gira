"""
Gathers and aggregates statistics on frequency of damage of transmission lines
"""

import os

TRANSMISSION_OUT = [
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
            STORM_IMPACT_STATISTICS_DIR, f"combined_storm_statistics_{config['central_threshold']}.csv"
        ),
    params:
        output_dir=config["output_dir"],
        aggregate_level=config["aggregate_level"],
        region_eval=STORM_BASINS,
        sample_eval=SAMPLES,
        nh_eval=STORMS,
        central_threshold=config["central_threshold"],
    output:
        TRANSMISSION_OUT,
    script:
        "./transmission_aggregate.py"
