# Open Global Infrastructure Risk/Resilience Analysis

[![mdBook Documentation](https://github.com/nismod/open-gira/actions/workflows/docs.yml/badge.svg?branch=main)](https://nismod.github.io/open-gira)
[![pyTest](https://github.com/nismod/open-gira/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/nismod/open-gira/actions/workflows/test.yml)

This open-source [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow will 
analyse physical climate risks to infrastructure networks using global open data. 

The related open-source Python library [snail](https://github.com/nismod/snail) provides
some of the core functionality.

> Work in Progress
>
> Goals: 
> - automated pipeline for reproducible analysis anywhere in the world
> - maps per-country and of larger areas
> - charts/stats of exposure per admin region, per hazard type, scenario, epoch
> - consider transport, electricity, water, communications systems
> - consider river flooding, storm surge coastal flooding, tropical cyclones
> - estimate direct damages to physical networks
> - estimate indirect effects of disruption - people affected, economic activity disrupted
>
> Non-goals:
> - will not build on closed data sources, which may be appropriate for other projects or use-cases
> - long-term planning or detailed operational simulation

## Installation

### conda

This repository comes with a `environment.yml` file describing the conda and pip packages required to run `open-gira`.

Create the `open-gira` conda environment:
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
python -m pytest tests
```

## Downloading datasets

New users can follow through the remaining steps in this guide using the Tanzania OpenStreetMap data, available from 
[https://download.geofabrik.de/africa/tanzania.html](https://download.geofabrik.de/africa/tanzania.html) (~500MB).
This file should be placed in ./data.

There should also be a hazard file in the location specified by `hazard_csv` in `config/config.yml`.
For an initial run, users can copy the files in `./tests/test_aqueduct_data` to the `<hazard_csv>` location 
(by default `./data/aqueduct`).

## Running the pipeline

The snakemake configuration details are in `config/config.yml`. 
You can edit this to set the target OSM
dataset, number of slices and hazard data location. See
[config/README.md](https://github.com/nismod/open-gira/blob/main/config/README.md)
for details on the configuration variables.
For new users, the default values should suffice.

The second step is to create a configuration file for `osmium
extract`, describing the bounding box of each OSM dataset slice.  This
configuration file is expected to be found next to the OSM data at the
location specified by the `data_dir` configuration variable. Its
filename should be of the form `<dataset>-extracts.json` where
`<dataset>` is specified by the `dataset` configuration variable. See
[Creating multiple extracts in one
go](https://osmcode.org/osmium-tool/manual.html#creating-geographic-extracts)
on the `osmium-tool` docs for more details on how to write config
files for `osmium extract`. A common task is to slice the OSM dataset
into areas of equal height and width, see [Automatically generating
the osmium extract confgiuration
file.](https://github.com/nismod/open-gira/tree/update_readme#step-by-step-description-of-the-pipeline)
New users should follow those steps to create `./tanzania-latests.json`.

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
## Automatically generating the `osmium extract` configuration file

A common task is to slice the OSM dataset into areas of equal height
and width. Script `prepare-extracts.py` automates this
process, given a JSON file describing the original dataset. 

Say that you want to slice `tanzania-latest.osm.pbf` into 6 slices of equal
height and equal width. First, write a `osmium extract` config file
describing the `tanzania-latest` as a single extract:

```json
// ./tanzania-latest.json
{
    "directory": "./results/slices",
    "extracts": [
        {
            "bbox": [
                -1.23,
                51.78,
                -1.175,
                51.805
            ],
            "output": "tanzania-latest.osm.pbf"
        }
	]
}
```

Next, use `prepare-extracts.py` to generate the `osmium extract`
configuration file for the 9 slices. For instance:

```
python prepare-extracts.py tanzania-latest.json 3
```

This generates a file `./data/tanzania-latest-extracts.json` describing
the 9 slices to be created by `osmium extract`.

## Step-by-step description of the pipeline

The pipeline starts from a OpenStreetMap dataset (_e.g._
`europe-latest`) and produces network/flood hazard intersection data,
associating road splits to corresponding flood levels.

The pipeline consists in the following steps:

1. The initial OSM dataset is filtered, keeping only relevant tags for road links
   (using `osmium tags-filter`. This results in a smaller file
   `<output_dir>/<dataset>.highway-core.osm.pbf`.
2. The OSM dataset's headers are examined for a `bbox` property and that is used
   to determine the bounding box for the whole area (`<output_dir>/json/<dataset>.json`).
3. The OSM dataset bounding box is sliced into a grid of smaller bounding boxes
   according to the `slice_count` config option (`<output_dir>/json/<dataset>-extracts.geojson`).
4. The filtered OSM file is sliced into areas of equal size using the bounding 
   box grid (`<output_dir>/slices/<dataset>-slice<N>.osm.pbf`).
5. Each filtered OSM dataset slice is then converted to the GeoParquet data format,
   resulting in `<output_dir>/geoparquet/<dataset>-slice<N>.highway-core.geoparquet`.
6. Each geoparquet slice is intersected against flood level data from the
   aqueduct dataset. The aqueduct dataset itself consists of a collection of
   raster data files. The network/hazard intersection results in data
   `<output_dir>/splits/<dataset>-slice<N>.highway-core.splits.geoparquet` describing
   roads split according to the raster grid and associated flood level values.
   A corresponding `parquet` files (without geometries) is also created.
7. Split data (one file per slice, see step 1) is then joined into a unique
   dataset describing splits and associated flood level values for the whole
   original OSM dataset. This results in
   `<output_dir>/<dataset>.highway-core_aqueduct_river_splits.geoparquet`.

### Keeping things tidy

You can remove all intermediate data and output files by running

```
snakemake clean
```

Note that this will *not* remove the final data files
`<output_dir>/<dataset>.highway-core.splits.[geoparquet, parquet]`.

Snakemake has utilities to improve the workflow code quality:
- `snakemake --lint` suggests improvements and fixes for common problems
- `snakefmt .` reformats files according to a code style guide, similar to `black` for Python code.

