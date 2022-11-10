# Process data

The processing of the [input data](../download.md) is split into the following steps
1. [Split world (world_splitter)](worldsplit.md)
2. [Document missing country data (process_exclude_countries)](excludecountries.md)
3. [Process powerplants (process_powerplants)](powerplants.md)
4. [Process targets (process_target_box)](targets.md)
5. [Process gridfinder (process_gridfinder)](gridfinder.md)
6. [Process network (process_network)](network.md)
7. [Process network connections (process_connector)](connector.md)

Here, the order is important to the processing and the snakemake workflow can be visualised as follows.

![`process_all` rule workflow with one box_id example](../img/dag_process_all.png)

Note that this includes the [download workflow](../download.md) process too. We also note that for this
visualisation only one box_id wildcard has been specified (box_1103) to prevent visual cluttering. Whilst this is
possible, the input usually consists of many more box_id wildcards which will be processed in parallel. This of course
depends on the config.yaml inputs.

These steps can be executed in one by entering `snakemake -s workflow/Snakefile process_all -c{core_numbers}`. Note
that if the `download_all` rule has not already run, then it will be run through this command automatically.
