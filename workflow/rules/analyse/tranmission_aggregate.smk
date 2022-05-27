"""Gathers and aggregates statistics on frequency of damage of transmission lines

"""
import os

transmission_out = [os.path.join(output_dir, "power_output", "statistics", "aggregate", "transmission_line_frequency_hit.gpkg"), os.path.join(output_dir, "power_output", "statistics", "aggregate", "transmission_line_reconstruction_costs.gpkg")]
rule analyse_transmission:
    input:
        os.path.join(stat_path, 'combined_storm_statistics.csv')
    params:
        output_dir = config['output_dir'],
        reconstruction_cost = config['reconstruction_cost'],
        aggregate_level = snakemake.params['aggregate_level'],
        region_eval = REGIONS,
        sample_eval = SAMPLES,
        nh_eval = config['specific_storm_analysis']
    output:
        transmission_out
    script:
        os.path.join('..', '..', 'scripts', 'analyse' ,'tranmission_aggregate.py')