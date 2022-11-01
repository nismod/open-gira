# Intersect infrastructure with storms

The intersection of the [storm data](../download/storm-ibtracs.md) with the
[processed infrastructure](../process.md) is ordered into the following steps
1. [Region boxes (`intersect_region_boxes`)](boxes.md)
2. [Generate region units (`intersect_unit_generator`)](grid.md)
3. [Analyse storm winds (`intersect_winds_indiv`)](windspeed.md)
4. [GDP losses (`intersect_damages`)](gdploss.md)

Here, the order is relevant. First the boxes are assigned their respective region (EP, NA, NI,
SI, SP, WP) in the [region boxing workflow](boxes.md). The next step ([unit
generation](grid.md)) generates 'units' in which the wind speeds are assumed to be the same
throughout. This is suitable for small units. Units containing no power network infrastructure
(i.e. no [network](../process/network.md)) are discarded.

A wind model is then applied to the [storm data](../download/storm-ibtracs.md) to
calculate the wind speed parameters at each unit ([wind analysis](windspeed.md)). Given a
failure criteria, infrastructure components are documented to fail in the [infrastructure
intersection workflow](gdploss.md).

The intersection (including previous workflow dependencies) can be visualised as follows.

![`intersect_all` rule workflow with one box_id, region and sample example](../../img/dag_intersect_all.png)

We use one box_id (box_1103), one region (NA) and one sample (0) for the wildcard parameters, however, in practice,
there are frequently many more. This is dependent on the config.yaml inputs.
The `intersect_all` rule has arrow dependencies other than `merge_overview_all_stats` due to the implementation of
snakemake checkpoints. This allows the correct scripts to be rerun in the event of incomplete previous runs (see
[wind analysis](windspeed.md) for further details).


These steps can be executed in one by entering

```
snakemake intersect_all -c{core_numbers}
```

Note that if the `process_all` and/or `download_all` rules have not already
run, then they will be run through this command automatically.

This process is used to bring some parameters of the storm intersection data together. It's main
purpose is actually such that snakemake is able to correctly use the checkpoint rule (as explained
in the [storm wind analysis](windspeed.md)). For each region and each sample, the
parameters described in the [damage calculation workflow](gdploss.md) are collected in
a region-sample csv file. Then each region-sample csv file is collected into a combined csv file `open-gira\results\power_output\statistics\combined_storm_statistics.csv`.
