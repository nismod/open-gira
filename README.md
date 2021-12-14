# OpenStreetMap / Aqueduct Flood exposure

A snakemake workflow based on `nismod/snail` to demonstrate infrastructure network
intersections with hazards as an exposure calculation.

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
dataset, number of slices and hazard data. See
[config/README.md](https://github.com/nismod/open-gira/blob/main/config/README.md)
for details on the configuration variables.

The second step is to create a configuration file for `osmium
extract`, describing the bounding box of each OSM dataset slice.  This
configuration file is expected to be found next to the OSM data at the
locatin specified by the `data_dir` configuration variable, with a
filename of the form `<dataset>-extracts.json` where `<dataset>` is
specified by the `dataset` configuration variable. See [Creating
multiple extracts in one
go](https://osmcode.org/osmium-tool/manual.html#creating-geographic-extracts)
on the `osmium-tool` docs for more details on how to write config
files for `osmium extract`. A common task is to slice the OSM dataset
into areas of equal height and width, see Automatically generating the
osmium extract confgiuratin file.

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
and width. Script `scripts/prepare-extracts.py` automates this
process, given a JSON file describing the original dataset. For
instance, say you want to slice `europe-latest.osm.pbf` into 6 slices
of equal height and equal width. First, write a `osmium extract`
config file describing the `europe-latest` as a single extract:

```json
{
    "directory": "./data",
    "extracts": [
        {
            "bbox": [
                -1.23,
                51.78,
                -1.175,
                51.805
            ],
		    "output": "europe-latest.osm.pbf"
        },
	]
}
```

Next, use `prepare-extracts.py` to generate the `osmium extract`
configuration file for the 6 slices. For instance:

```
python scripts/prepare-extracts.py europe-latest.json --width 3 --height 2
```

This generates a file `./data/europe-latest-extracts.json` describing
the 6 slices to be created by `osmium extract`.

## Step-by-step description of the pipeline

The pipeline starts from a OpenStreetMap dataset (_e.g._
`europe-latest`) and produces network/flood hazard intersection data,
associating road splits to corresponding flood levels.

The pipeline consists in the following steps:

1. The initial OSM dataset is sliced into areas of equal size
   (`results/slices/<dataset>-slice<N>.osm.pbf`).
2. Filters down each OSM data slice keeping only relevant tags for road links
   (using `osmium tags-filter`. This results in files
   `results/filtered/<dataset>-slice<N>.highway-core.osm.pbf`.
3. Each filtered OSM dataset is then converted to the GeoParquet data format,
   resulting in `results/geoparquet-slice<N>.highway-core.geoparquet`.
4. Each geoparquet slice is intersected against flood level data from the
   aqueduct dataset. The aqueduct dataset itself consists of a collection of
   raster data files. The network/hazard intersection results in data
   `results/splits/<dataset>-slice<N>.highway-core.splits.geoparquet` describing
   roads split according to the raster grid and associated flood level values.
   A corresponding `parquet` files (without geometries) is also created.
5. Split data (one file per slice, see step 1) is then joined into a unique
   dataset describing splits and associated flood level values for the whole
   original OSM dataset. This results in
   `results/<dataset>.highway-core.splits.geoparquet`.

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

