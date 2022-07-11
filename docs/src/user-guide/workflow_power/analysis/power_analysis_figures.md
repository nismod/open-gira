# Figure generation

In order to run the figure generation, all models (constant, CMCC-CM2-VHR4, CNRM-CM6-1-HR, EC-Earth3P-HR, HadGEM3-GC31-HM) must be run. 
To do so, once a model has run, rename the `results/power_output` folder to `results/power_output-{model}` e.g. `power_output-EC-Earth3P-HR`. Then change
the model in the config file and run `analyse_all` again. Repeat until you have all the following.
 - power_output-CMCC-CM2-VHR4
 - power_output-CNRM-CM6-1-HR
 - power_output-constant
 - power_output-EC-Earth3P-HR
 - power_output-HadGEM3-GC31-HM

Then run the `figures_all` rule to generate the figures. The outputs will be found in `results/power_figures`. To reiterate, you can NOT run the `figures_all` rule from the start, rather, first `analyse_all`, then the process as detailed above and then `figures_all`.
