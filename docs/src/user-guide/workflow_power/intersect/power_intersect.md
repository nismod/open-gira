# Intersect infrastructure with storms


The intersection of the [storm data](../download/power_download_stormtracks.md) with the
[processed infrastructure](../process/power_process.md) is ordered into the following steps 
1. [Region boxes (intersect_region_boxes)](power_intersect_boxes.md)
2. [Generate region units (intersect_unit_generator)](power_intersect_grid.md)
3. [Analyse storm winds (intersect_winds_indiv)](power_intersect_windextracter.md)
4. [Intersect with infrastructure (intersect_damages)](power_intersect_gdploss.md)
5. [Overview (merge_overview_indiv_stats and merge_overview_all_stats)](power_intersect_overview.md)

Here, the order is relevant. First the boxes are assigned their respective region (EP, NA, NI, SI, SP, WP) in the
[region boxing workflow](power_intersect_boxes.md). The next step ([unit generation](power_intersect_grid.md)) generates 'units' in which the wind speeds are assumed
to be the same throughout. This is suitable for small units. Units containing no power network infrastructure (i.e. no [network](../process/power_process_network.md)) are discarded. A wind
model is then applied to the [storm data](../download/power_download_stormtracks.md) to calculate the wind speed
parameters at each unit ([wind analysis](power_intersect_windextracter.md)). Given a failure criteria, infrastructure
components are documented to fail in the [infrastructure intersection workflow](power_intersect_gdploss.md). The
[overview](power_intersect_overview.md) brings together the storm intersections for a very preliminary overview.


The intersection (including previous workflow dependencies) can be visualised as follows.

![`intersect_all` rule workflow with one box_id, region and sample example](../power_img/dag_intersect_all.png)


We use one box_id (box_1103), one region (NA) and one sample (0) for the wildcard parameters, however, in practice,
there are frequently many more. This is dependent on the config.yaml inputs.
The `intersect_all` rule has arrow dependencies other than `merge_overview_all_stats` due to the implementation of
snakemake checkpoints. This allows the correct scripts to be rerun in the event of incomplete previous runs (see
[wind analysis](power_intersect_windextracter.md) for further details).


These steps can be executed in one by entering `snakemake -s workflow/Snakefile intersect_all -c{core_numbers}`. Note
that if the `process_all` and/or `download_all` rules have not already run, then they will be run through this command
automatically.