# Intersect Overview

This process is used to bring some parameters of the storm intersection data together. It's main
purpose is actually such that snakemake is able to correctly use the checkpoint rule (as explained
in the [storm wind analysis](power_intersect_windextracter.md)). For each region and each sample, the
parameters described in the [damage calculation workflow](power_intersect_gdploss.md) are collected in
a region-sample csv file. Then each region-sample csv file is collected into a combined csv file `open-gira\results\power_output\statistics\combined_storm_statistics.csv`.


