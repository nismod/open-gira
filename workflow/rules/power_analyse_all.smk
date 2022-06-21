"""Perform analysis of below selected workflows


"""


rule analyse_all:
    input:
        aggregate_levels_out,
        empirical_distribution_out,
        country_matrix_output,
        transmission_out,
        percentile_out,
        out_connector,
        completed_files  # to ensure storms are all run
