# Process data for power analysis

The processing of the [input data](../download/power_download.md) is split into the following steps
1. [Split world (world_splitter)](power_process_worldsplit.md)
2. [Document missing country data (process_exclude_countries)](power_process_excludecountries.md)
3. [Process powerplants (process_powerplants)](power_process_powerplants.md)
4. [Process targets (process_target_box)](power_process_targets.md)
5. [Process gridfinder (process_gridfinder)](power_process_gridfinder.md)
6. [Process network (process_network)](power_process_network.md)
7. [Process network connections (process_connector)](power_process_connector.md)

Here, the order is important to the processing and the snakemake workflow can be visualised as follows.

![`process_all` rule workflow with one box_id example](../power_img/dag_process_all.png)

Note that this includes the [download workflow](../download/power_download.md) process too. We also note that for this visualisation only one box_id wildcard has been specified (box_1103) to prevent visual cluttering. Whilst this is possible, the input usually consists of many more box_id wildcards which will be processed in parallel. This of course depends on the config.yaml inputs (see more details [here](../../workflow_power/workflow_power.md))