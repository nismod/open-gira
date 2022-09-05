"""
Perform analysis of below selected workflows
"""


rule analyse_all:
    input:
        AGGREGATE_LEVELS_OUT,
        EMPIRICAL_DISTRIBUTION_OUT,
        COUNTRY_MATRIX_OUTPUT,
        TRANSMISSION_OUT,
        PERCENTILE_OUT,
        CONNECTOR_OUT,
        COMPLETION_FLAG_FILES,  # to ensure storms are all run
