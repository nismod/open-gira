# Configuration

The pipeline is configured from `config/config.yaml` file. 
We provide a basic config file `config/config.yml` with explanatory comments that you can modify.

The configuration file is meant to specify the location of input data
and outputs, as well as initial OSM dataset and its slicing.

- `data_dir`: Relative or absolute path to directory containing initial OSM dataset.
- `hazard_data_dir`: Relative or absolute path to directory containing aqueduct data files.
- `output_dir`: Relative or absolute path where the output should be placed. Will be created if it does not exist.

- `hazard_csv`: Name of CSV file describing aqueduct dataset. Example:

  |key|climate\_scenario|model|year|return\_period|filename|
  |---|----------------|-----|----|-------------|--------|
  |0|inunriver\_rcp8p5\_00IPSL-CM5A-LR\_2080\_rp00050|rcp8p5|IPSL-CM5A-LR|2080|50|inunriver\_rcp8p5\_00IPSL-CM5A-LR\_2080\_rp00050.tif|
  |1|inunriver\_rcp4p5\_00000NorESM1-M\_2030\_rp00005|rcp4p5|NorESM1-M|2030|5|inunriver\_rcp4p5\_00000NorESM1-M\_2030\_rp00005.tif|
  |...|...|...|...|...|...|

  Filename should be relative to `<hazard_data_dir>`.
- `dataset`: Name of initial OpenStreetMap dataset, _e.g._ `tanzania-latest`.
- `edge_attrs`: Edge attributes list for network/hazard intersection.
- `osmium_tags_filters_file`: File containing the OSM attributes to filter the input data

Modifying the configuration file will *not* trigger a re-run of the pipeline by
snakemake. If you wish to rerun the whole pipeline after altering the
configuration, use 

```
snakemake -R
```

Re-running the whole pipeline from the start might not be necessary. For
instance if you modify `<output_dir>`, only the last 2 pipeline steps will be
altered. In this case, you can ask snakemake to (re)start from the first
affected rule (see `Snakefile`), and it will figure out what must be done to
complete the pipeline. In this case: 

```
snakemake -R network_hazard_intersection
```

which will re-run `network_hazard_intersection` and `join_data` rule,
in this order.
