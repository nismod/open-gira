"""
Plots the empirical distribution of storms for simple statistics.
"""

import os


EMPIRICAL_DISTRIBUTION_OUT = [os.path.join(STORM_IMPACT_STATISTICS_DIR, 'empirical', f'empirical_{metric}.png') for metric in STORM_ANALYSIS_METRICS]


rule analyse_empirical_distribution:
    input:
        STORM_STATS_BY_THRESHOLD
    params:
        output_dir = config['output_dir'],
        metrics = STORM_ANALYSIS_METRICS,
        region_eval = STORM_BASINS,
        sample_eval = SAMPLES,
        nh_eval = STORMS,
        central_threshold = config['central_threshold'],
        minimum_threshold = config['minimum_threshold'],
        maximum_threshold = config['maximum_threshold']
    output:
        EMPIRICAL_DISTRIBUTION_OUT
    script:
        os.path.join('..', '..', 'scripts', 'analyse' , 'storm_distribution_empirical.py')
