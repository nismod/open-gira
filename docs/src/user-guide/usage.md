# Usage

`open-gira` is comprised of a set of `snakemake` rules which call scripts and
library code to request data, process it and produce results.

## Snakemake

The key idea of `snakemake` is similar to `make` in that the workflow is
determined from the end (the files users want) to the beginning (the files
users have, if any) by applying general rules with pattern matching on file and
folder names.

A example invocation looks like:

```bash
snakemake --cores 2 -- results/wales-latest_filter-road-primary/edges.gpq
```

Here, we ask `snakemake` to use up to 2 CPUs to produce a target file, in this
case, the edges of the Welsh road network. `snakemake` pattern matches
`wales-latest` as the OSM dataset name and `filter-road` as the network type we
want to filter for.

To check what work we're going to request before commencing, use the `-n` flag:

```bash
snakemake -n --cores 2 -- results/wales-latest_filter-road-primary/edges.gpq
```

This will explain which rules will be required to run to produce the target
file. It may be helpful to [visualise](https://snakemake.readthedocs.io/en/stable/executing/cli.html#visualization)
which rules are expected to run, too.

## Configuration

The snakemake configuration details are in `config/config.yml`. You can edit
this to set the target OSM infrastructure datasets, number of spatial slices, and
hazard datasets. See
[config/README.md](https://github.com/nismod/open-gira/blob/main/config/README.md)
and the documentation for each workflow for more details on the configuration
variables.
