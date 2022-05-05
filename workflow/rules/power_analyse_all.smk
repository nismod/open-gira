"""Perform analysis of below selected workflows


"""


rule analyse_all:
    input:
        aggregate_levels_out,
        empirical_distribution_out,
        country_matrix_output
