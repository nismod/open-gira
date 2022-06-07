"""Plots the empirical distribution of storms for simple statistics.

"""
import os
empirical_distribution_out = [os.path.join(stat_path, 'empirical', f'empirical_{metric}.png') for metric in metrics]

rule analyse_empirical_distribution:
    input:
        stat_csv
    params:
        output_dir = config['output_dir'],
        metrics = metrics,
        region_eval = REGIONS,
        sample_eval = SAMPLES,
        nh_eval = STORMS,
        central_threshold = config['central_threshold'],
        minimum_threshold = config['minimum_threshold'],
        maximum_threshold = config['maximum_threshold']
    output:
        empirical_distribution_out
    script:
        os.path.join('..', '..', 'scripts', 'analyse' ,'storm_distribution_empirical.py')