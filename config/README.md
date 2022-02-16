# Configuration

The pipeline is configured from `config/config.yaml` file. 
We provide a basic config file `config/config.yml` with explanatory comments that you can modify.

The configuration file is meant to specify the location of input data
and outputs, as well as some runtime settings.

- `output_dir`: Relative or absolute path where the output should be placed. Will be created if it does not exist.
- `hazard_datasets`: Named list of file paths to `.txt` files containing a list of hazard files.
Files can be specified (both the `.txt` files and the `.tif` files they point to) either as filenames
relative to the project root, as absolute file paths, or as remote resources (`http://` or `https://`). 
The names in the list should not include `_` or `/` characters.
Remote resources will be fetched with the `wget` utility.
- `infrastructure_datasets`: Named list of file paths to `.osm.pbf` files to use as datasets. 
These can be local files, specified with absolute file paths or file paths relative to the project root,
or they can be remote files to fetch with `wget`.

- `slice_count`: Number of slices to take for each infrastructure dataset.
More slices allows for greater parallelization, but will also duplicate sections of roads that cross
slice boundaries, so too high a number can lead to redundancy.
- `edge_attrs`: Edge attributes list for network/hazard intersection.
- `osmium_tags_filters_file`: File containing the OSM attributes to filter the input data

Modifying the configuration file will *not* trigger a re-run of the pipeline by
snakemake. If you wish to rerun the whole pipeline after altering the
configuration, use 

```
snakemake --cores all --forceall
```

Re-running the whole pipeline from the start might not be necessary. For
instance if you modify `<output_dir>`, only the last few pipeline steps will be
altered. In this case, you can ask snakemake to (re)start from the first
affected rule (see `Snakefile`), and it will figure out what must be done to
complete the pipeline. In this case: 

```
snakemake --cores all -R intersection
```

which will re-run `intersection` and `join_data` rules,
in this order.
