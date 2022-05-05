""" Contains the common data for the analysis rules"""

import os
stat_path = os.path.join(config['output_dir'], 'power_output', 'statistics')
metrics = ['GDP losses', 'targets with no power (f=0)', 'population affected', 'population with no power (f=0)', 'effective population affected']
metrics_target = ['population', 'mw_loss_storm', 'f_value', 'gdp_damage']
increased_severity_sort_bool = str(config["increased_severity_sort"])[0]