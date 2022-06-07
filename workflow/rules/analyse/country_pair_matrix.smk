"""Plots the empirical storm relationship matrix between (two) countries and conditional probability relationship

"""
import os

title_unint = "Empirical fraction of a storm hitting both countries given it hits one of them"
title_condprob = "Given country A is hit, what is the likelihood also B is hit"

country_matrix_output = [os.path.join(stat_path, 'empirical', title_unint+".png"), os.path.join(stat_path, 'empirical', title_condprob+".png")]


rule analyse_country_matrix:
    input:
        os.path.join(stat_path, f"combined_storm_statistics_{config['central_threshold']}.csv")
    params:
        output_dir = config['output_dir'],
        central_threshold = config['central_threshold']
    output:
        country_matrix_output
    script:
        os.path.join('..', '..', 'scripts', 'analyse' ,'storm_distribution_empirical_country_matrix.py')