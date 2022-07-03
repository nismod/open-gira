"""for all figures

"""
import os


rule figures_all:
    input:
        initial_checks,
        out_diff_EACA_file,
        out_diff_EAD_agg,
        out_diff_EAD_agg_norm,
        out_diff_EAD_indiv,
        out_cc_file_diff,
        out_cc_file,
        out_cdf_EAD,
        out_cdf_EACA,
        out_agg_EACA_plot,
        out_agg_EACA_plot_perc,
        out_diff_EAD_indiv_plot,
        out_current_EAD_indiv_plot,
        #out_diff_EAD_plot,
        #out_diff_EAD_plot_norm
