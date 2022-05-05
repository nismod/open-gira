"""Plots the empirical distribution of storms for simple statistics.

"""
import os
empirical_distribution_out = [os.path.join(stat_path, 'empirical', f'empirical_{metric}.png') for metric in metrics]

rule analyse_empirical_distribution:
    input:
        os.path.join(stat_path, 'combined_storm_statistics.csv')
    params:
        output_dir = config['output_dir'],
        metrics = metrics
    output:
        empirical_distribution_out
    script:
        os.path.join('..', '..', 'scripts', 'analyse' ,'storm_distribution_empirical.py')