"""Aggregates  # TODO



constant (no climate effects), CMCC-CM2-VHR4, CNRM-CM6-1-HR, EC-Earth3P-HR, HadGEM3-GC31-HM


"""
import os


## aggregated (absolute) ##
metric_EAD_agg = 'reconstruction_cost_annual_expected'
merge_key_EAD_agg = 'code'

## aggregated but normalised ##
metric_EAD_agg_norm = 'reconstruction_cost_annual_expected_fraction_normalised'
merge_key_EAD_agg_norm = 'code'

## individual assets ##
metric_EAD_indiv = 'reconstruction_cost_annual_expected'
merge_key_EAD_indiv = 'link'



### outputs ###
out_diff_EAD_agg = os.path.join(config['output_dir'], 'power_figures', 'intermediate_files', f'difference_{merge_key_EAD_agg}_{metric_EAD_agg}.gpkg')
out_diff_EAD_agg_norm = os.path.join(config['output_dir'], 'power_figures', 'intermediate_files', f'difference_{merge_key_EAD_agg_norm}_{metric_EAD_agg_norm}.gpkg')
out_diff_EAD_indiv = os.path.join(config['output_dir'], 'power_figures', 'intermediate_files', f'difference_{merge_key_EAD_indiv}_{metric_EAD_indiv}.gpkg')




## aggregated (absolute) ##
rule fig_aggregate_EAD_agg:
    input:
        [os.path.join(config['output_dir'], f'power_output-{model}', 'statistics', 'aggregate', "transmission_line_reconstruction_costs.gpkg") for model in models_future]
    params:
        output_dir = config['output_dir'],
        metric = metric_EAD_agg,
        merge_key = merge_key_EAD_agg
    output:
        out_agg_EAD_agg = os.path.join(config['output_dir'], 'power_figures', 'intermediate_files', f'mean_{merge_key_EAD_agg}_{metric_EAD_agg}.gpkg')
    script:
        os.path.join("..", "..", "..", 'scripts', 'analyse', 'figures', 'mean_agg.py')


rule fig_diff_EAD_agg:
    input:
        [rules.fig_aggregate_EAD_agg.output.out_agg_EAD_agg, os.path.join(config['output_dir'], f'power_output-constant', 'statistics', 'aggregate', "transmission_line_reconstruction_costs.gpkg")]  # first file is future
    params:
        metric = metric_EAD_agg,
        merge_key = merge_key_EAD_agg
    output:
        out_diff_EAD_agg
    script:
        os.path.join("..", "..", "..", 'scripts', 'analyse', 'figures', 'diff_agg.py')



## aggregated but normalised ##
rule fig_aggregate_EAD_agg_norm:
    input:
        [os.path.join(config['output_dir'], f'power_output-{model}', 'statistics', 'aggregate', "transmission_line_reconstruction_costs.gpkg") for model in models_future]
    params:
        output_dir = config['output_dir'],
        metric = metric_EAD_agg_norm,
        merge_key = merge_key_EAD_agg_norm
    output:
        out_agg_EAD_agg_norm = os.path.join(config['output_dir'], 'power_figures', 'intermediate_files', f'mean_{merge_key_EAD_agg_norm}_{metric_EAD_agg_norm}.gpkg')
    script:
        os.path.join("..", "..", "..", 'scripts', 'analyse', 'figures', 'mean_agg.py')


rule fig_diff_EAD_agg_norm:
    input:
        [rules.fig_aggregate_EAD_agg_norm.output.out_agg_EAD_agg_norm, os.path.join(config['output_dir'], f'power_output-constant', 'statistics', 'aggregate', "transmission_line_reconstruction_costs.gpkg")]  # first file is future
    params:
        metric = metric_EAD_agg_norm,
        merge_key = merge_key_EAD_agg_norm
    output:
        out_diff_EAD_agg_norm
    script:
        os.path.join("..", "..", "..", 'scripts', 'analyse', 'figures', 'diff_agg.py')




## individual assets ##
rule fig_aggregate_EAD_indiv:
    input:
        [os.path.join(config['output_dir'], f'power_output-{model}', 'statistics', 'aggregate', "transmission_line_frequency_hit.gpkg") for model in models_future]
    params:
        output_dir = config['output_dir'],
        metric = metric_EAD_indiv,
        merge_key = merge_key_EAD_indiv
    output:
        out_agg_EAD_indiv = os.path.join(config['output_dir'], 'power_figures', 'intermediate_files', f'mean_{merge_key_EAD_indiv}_{metric_EAD_indiv}.gpkg')
    script:
        os.path.join("..", "..", "..", 'scripts', 'analyse', 'figures', 'mean_agg.py')


rule fig_diff_EAD_indiv:
    input:
        [rules.fig_aggregate_EAD_indiv.output.out_agg_EAD_indiv, os.path.join(config['output_dir'], f'power_output-constant', 'statistics', 'aggregate', "transmission_line_frequency_hit.gpkg")]  # first file is future
    params:
        metric = metric_EAD_indiv,
        merge_key = merge_key_EAD_indiv
    output:
        out_diff_EAD_indiv
    script:
        os.path.join("..", "..", "..", 'scripts', 'analyse', 'figures', 'diff_agg.py')

