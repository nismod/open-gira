# Intersect Overview

This process is used to bring some parameters of the storm intersection data together. It's main
purpose is actually such that snakemake is able to correctly use the checkpoint rule (as explained
in the [storm wind analysis](power_intersect_windextracter.md)). For each region and each sample, the
parameters described in the [damage calculation workflow](power_intersect_gdploss.md) are collected in
a region-sample csv file. Then each region-sample csv file is collected into a combined csv file `open-gira\results\power_output\statistics\combined_storm_statistics.csv`.


### Manual run
It is possible to manually run this file by calling `python3 open-gira\workflow\scripts\intersect\intersect_overview_individual.py` first and then
`python3 open-gira\workflow\scripts\intersect\intersect_overview.py`. These files will manually scan for processed storms
and so can be prone to issues if files are misplaced or missing. The output of the latter command is `open-gira\results\power_output\statistics\combined_storm_statistics_manual.csv`.