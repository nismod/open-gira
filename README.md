# OpenStreetMap / Aqueduct Flood exposure

Use `nismod/snail` to demonstrate infrastructure network intersections with hazards as an
exposure calculation.

Goals: 
- automated pipeline for reproducible exposure analysis anywhere in the world.
- map of a large area showing exposure at different return periods
- charts/stats of exposure per admin region (by road/rail) per hazard type, scenario, epoch

## Installation

### conda

This repository comes with a `environment.yml` file describing the conda and pip packages required to run `open-gira`.

First, create the `open-gira` conda environment:
```
conda env create -f workflow/envs/environment.yml
```
and activate it
```
conda activate open-gira
```

### pip

Install python requirements as listed in `requirements.txt` - for
example using a venv:

```
python3 -m venv ./venv
. venv/bin/activate
pip install -r requirements.txt
```

### Install osmium-tool

Install [`osmium-tool`](https://osmcode.org/osmium-tool/manual.html) according to the
instructions there. Tests run with versions:
- osmium-tool v1.13.2
- libosmium v2.17.1

## Running tests

Workflow steps are tested using a small sample dataset. Run:

```
python -m pytest .tests
```


## Running the pipeline

Start by making your own copy of `config/config_template.yml` named
`config/config.yml`. You can edit the latter to set the target OSM
dataset, number of slices and hazard data. See **??** for details on
the configuration.

You can then run the exposure analysis pipeline automatically using
snakemake, like so

```
snakemake --cores 8
```

Individual configuration parameters can be overriden from the command
line, for instance

```
snakemake --cores 1 --config dataset=tanzania-latest data_dir=osm-data
```

It is often useful to maintain several configuration files. You can
specify a configuration to be used in place of the default
`config/config.yml` like so

```
snakemake --cores 8 --configfile config/my_other_config.yml
```


The pipeline consists in the following steps:

1. Slices the initial OSM dataset into areas of equal size
   (`<data_dir>/<dataset>-slice<N>.osm.pbf`).
2. Filters down each OSM data slice keeping only relevant tags for road links
   (using `osmium tags-filter`. This results in files
   `<data_dir>/<dataset>-slice<N>.highway-core.osm.pbf`.
3. Each filtered OSM dataset is then converted to the GeoParquet data format,
   resulting in `<data_dir>/<dataset>-slice<N>.highway-core.geoparquet`.
4. Each geoparquet slice is intersected against flood level data from the
   aqueduct dataset. The aqueduct dataset itself consists of a collection of
   raster data files. The network/hazard intersection results in data
   `<output_dir>/<dataset>-slice<N>.highway-core.splits.geoparquet` describing
   roads split according to the raster grid and associated flood level values.
   A corresponding `parquet` files (without geometries) is also created.
5. Split data (one file per slice, see step 1) is then joined into a unique
   dataset describing splits and associated flood level values for the whole
   original OSM dataset. This results in
   `<dataset>.highway-core.splits.geoparquet`.

### Configuration

The pipeline can be configured providing a `config.yaml` file. It is meant to
specify the location of input data and outputs, as well as initial OSM dataset
and its slicing.

- `data_dir`: Relative or absolute path to directory containing initial OSM dataset.
- `aqueduct_dir`: Relative or absolute path to directory containing aqueduct data files.
- `output_dir`: Relative or absolute path to directory containing output files (split data).
- `datafiles_list`: Name of CSV file describing aqueduct dataset. Example:

  |key|climate\_scenario|model|year|return\_period|filename|
  |---|----------------|-----|----|-------------|--------|
  |0|inunriver\_rcp8p5\_00IPSL-CM5A-LR\_2080\_rp00050|rcp8p5|IPSL-CM5A-LR|2080|50|inunriver\_rcp8p5\_00IPSL-CM5A-LR\_2080\_rp00050.tif|
  |1|inunriver\_rcp4p5\_00000NorESM1-M\_2030\_rp00005|rcp4p5|NorESM1-M|2030|5|inunriver\_rcp4p5\_00000NorESM1-M\_2030\_rp00005.tif|
  |...|...|...|...|...|...|

  Filename should be relative to `<aqueduct_dir>`.
- `dataset`: Name of initial OpenStreetMap dataset, _e.g._ `spain-latest`.
- `ratio`: Ratio for slicing original dataset. A ratio of 3 will
  result in 9 slices of equal area.

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

### Keeping things tidy

You can remove all intermediate data and output files by running

```
snakemake clean
```

Note that this will *not* remove the final data files
`<output_dir>/<dataset>.highway-core.splits.[geoparquet, parquet]`.

## Data

### Gridfinder

```
pip install zenodo_get mkdir -p data/gridfinder cd data/gridfinder
zenodo_get -d 10.5281/zenodo.3628142
```

### WRI Aqueduct

Use `nismod/aqueduct`.

TODO update to provide snail raster dataset metadata.

### OpenStreetMap

Testing with downloads from http://download.geofabrik.de/asia.html - Laos and Bangladesh.

TODO consider planet and diffs.

