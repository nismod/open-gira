"""Aggregates  # TODO



constant (no climate effects), CMCC-CM2-VHR4, CNRM-CM6-1-HR, EC-Earth3P-HR, HadGEM3-GC31-HM


"""
import os

increased_severity_sort_bool = str(config['increased_severity_sort'])[0]
assert increased_severity_sort_bool in ['T', 'F']

metric_EACA = 'effective_population_anually-expected'
merge_key_EACA = 'code'

out_agg_EACA_file = os.path.join(config['output_dir'], 'figures', 'intermediate_files', f'mean_{merge_key_EACA}_{metric_EACA}.gpkg')
out_diff_EACA_file = os.path.join(config['output_dir'], 'figures', 'intermediate_files', f'difference_{merge_key_EACA}_{metric_EACA}.gpkg')
print(os.path.join('workflow', 'scripts', 'analyse', 'figures', 'mean_agg.py'))
rule fig_aggregate_EACA:
    input:
        [os.path.join(config['output_dir'], f'power_output-{model}', 'statistics', 'aggregate', f"targets_geo_top{config['top_select']}{increased_severity_sort_bool}percent_aggregated_region.gpkg") for model in models_future]
    params:
        output_dir = config['output_dir'],
        metric = metric_EACA,
        merge_key = merge_key_EACA
    output:
        out_agg_EACA = out_agg_EACA_file
    script:
        os.path.join("..", "..", "..", "scripts", "analyse", "figures", "mean_agg.py")


rule fig_diff_EACA:
    input:
        [rules.fig_aggregate_EACA.output.out_agg_EACA, os.path.join(config['output_dir'], f'power_output-constant', 'statistics', 'aggregate', f"targets_geo_top{config['top_select']}{increased_severity_sort_bool}percent_aggregated_region.gpkg")]  # first file is future
    params:
        output_dir = config['output_dir'],
        metric = metric_EACA,
        merge_key = merge_key_EACA
    output:
        out_diff_EACA_file
    script:
        os.path.join("..", "..", "..", 'scripts', 'analyse', 'figures', 'diff_agg.py')
