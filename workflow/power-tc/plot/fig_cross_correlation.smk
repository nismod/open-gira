"""Plots the cross correlation JHR and future differences

"""
import os


in_cc = [
    os.path.join(
        config["output_dir"],
        f"power_output-{model}",
        "statistics",
        "empirical",
        "empirical_plotting_data",
        "country_matrix_both.csv",
    )
    for model in models_all
]

out_cc_file = [
    os.path.join(config["output_dir"], "power_figures", f"JHR_{name_cc_constant}.png"),
    os.path.join(config["output_dir"], "power_figures", f"JHR_{name_cc_future}.png"),
]
out_cc_file_diff = [
    os.path.join(
        config["output_dir"], "power_figures", f"JHR_{name_cc_future_diff}.png"
    ),
    os.path.join(
        config["output_dir"], "power_figures", f"JHR_{name_cc_future_perc_diff}.png"
    ),
]


rule fig_cross_correlation:
    input:
        in_cc,
    params:
        output_dir=config["output_dir"],
        remove_countries=remove_countries,
        name_cc_future=name_cc_future,
        name_cc_constant=name_cc_constant,
    output:
        out_cc=out_cc_file,
    script:
        "./cross_correlation.py"


rule fig_cross_correlation_diff:
    input:
        in_cc,
    params:
        output_dir=config["output_dir"],
        remove_countries=remove_countries,
        name_cc_future_diff=name_cc_future_diff,
        name_cc_future_perc_diff=name_cc_future_perc_diff,
    output:
        out_cc_file_diff,
    script:
        "./diff_cross_correlation.py"
