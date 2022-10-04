"""
Plots the empirical storm relationship matrix between (two) countries and conditional probability relationship
"""

import os


COUNTRY_MATRIX_OUTPUT = [
    os.path.join(STORM_IMPACT_STATISTICS_DIR, "empirical", "storms_hitting_A_and_B_given_either_hit" + ".png"),
    os.path.join(STORM_IMPACT_STATISTICS_DIR, "empirical", "likelihood_B_hit_given_A_hit" + ".png"),
]


rule analyse_country_matrix:
    conda: "../../../environment.yml"
    input:
        os.path.join(
            STORM_IMPACT_STATISTICS_DIR, f"combined_storm_statistics_{config['central_threshold']}.csv"
        ),
    params:
        output_dir=config["output_dir"],
        central_threshold=config["central_threshold"],
    output:
        COUNTRY_MATRIX_OUTPUT,
    script:
        os.path.join(
            "..",
            "..",
            "scripts",
            "analyse",
            "storm_distribution_empirical_country_matrix.py",
        )
