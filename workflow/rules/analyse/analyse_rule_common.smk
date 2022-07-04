""" Contains the common data for the analysis rules"""

import os

stat_path = os.path.join(config["output_dir"], "power_output", "statistics")
metrics = [
    "GDP losses",
    "targets with no power (f=0)",
    "population affected",
    "population with no power (f=0)",
    "effective population affected",
    "reconstruction cost",
]  # metrics general storm
metrics_target = [
    "population_without_power",
    "effective_population",
    "affected_population",
    "mw_loss_storm",
    "f_value",
    "gdp_damage",
]  # metrics targets
increased_severity_sort_bool = str(config["increased_severity_sort"])[0]
