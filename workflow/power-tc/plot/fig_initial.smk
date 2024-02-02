"""Checks this part of workflow can be proceeded with and provides common variables

Reminder: do not run figures_all until analyse_all has been completed!
"""


import os

# prelim variables
models_future = ["CMCC-CM2-VHR4", "CNRM-CM6-1-HR", "EC-Earth3P-HR", "HadGEM3-GC31-HM"]
models_all = ["constant"] + models_future  # constant MUST be first
remove_countries = ["USA", "VEN", "CYM", "VCT", "BHS", "ATG", "DMA", "LCA", "TTO"]
name_cc_constant = "constant climate"
name_cc_future = "future climate"
name_cc_future_diff = "Future minus current"
name_cc_future_perc_diff = "Future minus current (change in percent of current)"

master_linewidth = 0.55  # linewidth when lines in GeoDataFrame


initial_checks = os.path.join(
    config["output_dir"], "power_figures", "initial_checks.txt"
)


rule fig_checks:
    params:
        output_dir=config["output_dir"],
        models_all=models_all,
    output:
        initial_checks,
    script:
        "./fig_initial.py"
