"""Aggregates EAD and finds difference. For indivudal asset and aggregate

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
metric_EAD_indiv_perc = 'perc_'+metric_EAD_indiv
metric_EAD_indiv_perkm = metric_EAD_indiv+'_per_km'
merge_key_EAD_indiv = 'link'



### outputs ###
out_diff_EAD_agg = os.path.join(config['output_dir'], 'power_figures', 'intermediate_files', f'difference_{merge_key_EAD_agg}_{metric_EAD_agg}.gpkg')
out_diff_EAD_agg_norm = os.path.join(config['output_dir'], 'power_figures', 'intermediate_files', f'difference_{merge_key_EAD_agg_norm}_{metric_EAD_agg_norm}.gpkg')
out_diff_EAD_indiv = os.path.join(config['output_dir'], 'power_figures', 'intermediate_files', f'difference_{merge_key_EAD_indiv}_{metric_EAD_indiv}.gpkg')

out_diff_EAD_indiv_plot = os.path.join(config['output_dir'], 'power_figures', 'EAD_difference_perc_recon.png')
out_current_EAD_indiv_plot = os.path.join(config['output_dir'], 'power_figures', 'EAD_current_recon_per_km.png')
out_diff_EAD_plot = os.path.join(config['output_dir'], 'power_figures', 'EAD_difference_recon.png')
out_diff_EAD_plot_norm = os.path.join(config['output_dir'], 'power_figures', 'EAD_difference_recon_norm.png')

## aggregated (absolute) ##
rule fig_aggregate_EAD_agg:
    input:
        [os.path.join(config['output_dir'], f'power_output-{model}', 'statistics', 'aggregate', "transmission_line_reconstruction_costs.gpkg") for model in models_future],
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


# rule fig_diff_EAD_agg_plot:
#       """Not ideal plot as varies a lot by total km in region"""
#     input:
#         out_diff_EAD_agg
#     params:
#         output_dir = config['output_dir'],
#         metric = metric_EAD_agg,
#         vmax = 500000,
#         vmin = -500000,
#         cmap = 'RdBu_r',
#         linewidth = None
#     output:
#         out_diff_EAD_plot
#     script:
#         os.path.join("..", "..", "..", 'scripts', 'analyse', 'figures', 'plotter.py')




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

# rule fig_diff_EAD_agg_plot_norm:
#     input:
#         out_diff_EAD_agg
#     params:
#         output_dir = config['output_dir']
#         metric = metric_EAD_agg_norm,
#         vmax = 5000,
#         vmin = -5000,
#         cmap = 'RdBu_r',
#         linewidth = None
#     output:
#         out_diff_EAD_plot_norm
#     script:
#         os.path.join("..", "..", "..", 'scripts', 'analyse', 'figures', 'plotter.py')
#



## individual assets ##
rule fig_aggregate_EAD_indiv:
    input:
        in_agg_EAD_indiv = [os.path.join(config['output_dir'], f'power_output-{model}', 'statistics', 'aggregate', "transmission_line_frequency_hit.gpkg") for model in models_future]
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



rule fig_diff_EAD_indiv_plot:
    """Plots difference"""
    input:
        out_diff_EAD_indiv
    params:
        output_dir = config['output_dir'],
        metric = metric_EAD_indiv_perc,
        vmax = 100,
        vmin = -100,
        cmap = 'RdBu_r',
        linewidth = master_linewidth,
        legend_name = 'EAD change per km [% of constant climate]'
    output:
        out_diff_EAD_indiv_plot
    script:
        os.path.join("..", "..", "..", 'scripts', 'analyse', 'figures', 'plotter.py')

rule fig_current_EAD_indiv_plot:
    """Plots current"""
    input:
        os.path.join(config['output_dir'], f'power_output-constant', 'statistics', 'aggregate', "transmission_line_frequency_hit.gpkg")
    params:
        output_dir = config['output_dir'],
        metric = metric_EAD_indiv_perkm,
        vmax = 250,
        vmin = 0,
        cmap = 'Reds',
        linewidth = master_linewidth,
        legend_name = 'EAD per km [USD]'
    output:
        out_current_EAD_indiv_plot
    script:
        os.path.join("..", "..", "..", 'scripts', 'analyse', 'figures', 'plotter.py')


