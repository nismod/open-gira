# Open Global Infrastructure Risk/Resilience Analysis

[![mdBook Documentation](https://github.com/nismod/open-gira/actions/workflows/docs.yml/badge.svg?branch=main)](https://nismod.github.io/open-gira)
[![pyTest](https://github.com/nismod/open-gira/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/nismod/open-gira/actions/workflows/test.yml)
[![snakemake workflow](https://img.shields.io/badge/snakemake-open--gira-informational)](https://snakemake.github.io/snakemake-workflow-catalog/?usage=nismod/open-gira)

This open-source [snakemake](https://snakemake.readthedocs.io/en/stable/)
workflow can be used to analyse physical climate risks to infrastructure
networks using global open data.

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

### Conda packages

This repository comes with a `environment.yml` file describing the `conda` and
`PyPI` packages required to run `open-gira`. The `open-gira` developers
recommend using either [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html#micromamba)
or [mamba](https://mamba.readthedocs.io/en/latest/index.html) to install and
manage these `conda` packages.

#### Locally

Having installed one of the suggested package managers, to create the
`open-gira` conda environment:
```bash
micromamba create -f environment.yml -y
```

And to activate the environment:
```bash
micromamba activate open-gira
```

You are now ready to request result files, triggering analysis jobs in the
process.

#### Cluster

If installing on a cluster, you can work as above, or, create a seperate
orchestrating environment containing only snakemake, e.g.
```bash
micromamba create -n snakemake python=3.9 snakemake
```

In this context, `snakemake` itself can manage the other required dependencies,
creating other environments as necessary. To activate the orchestration
environment:
```bash
micromamba activate snakemake
```

To run the workflow on a [cluster](https://snakemake.readthedocs.io/en/stable/executing/cluster.html)
you will need to provide a [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles),
requesting targets as follows:
```bash
snakemake --profile <path_to_cluster_config> -- <target_file>
```

### exactextract

`exactextract` is used for zonal statistics in the tropical cyclones /
electricity grid analysis. It is not available via the `conda` package
management ecosystem and so must be installed separately. Please see
installation instructions [here](https://github.com/isciences/exactextract).

## Running tests

Workflow steps are tested using a small sample dataset. To run the tests:
```bash
python -m pytest tests
```

## Usage

`open-gira` is comprised of a set of `snakemake` rules which call scripts and
library code to request data, process it and produce results.

The key idea of `snakemake` is similar to `make` in that the workflow is
determined from the end (the files users want) to the beginning (the files
users have, if any) by applying general rules with pattern matching on file and
folder names.

A example invocation looks like:
```
snakemake --cores 2 -- results/wales-latest_filter-road/edges.geoparquet
```

Here, we ask `snakemake` to use up to 2 CPUs to produce a target file, in this
case, the edges of the Welsh road network. To check what work we're going to
request before commencing, use the `-n` flag:
```
snakemake -n --cores 2 -- results/wales-latest_filter-road/edges.geoparquet
```

This will explain which rules will be required to run to produce the target
file. It may be helpful to [visualise](https://snakemake.readthedocs.io/en/stable/executing/cli.html#visualization)
which rules are expected to run, too.

### Configuration

The snakemake configuration details are in `config/config.yml`. You can edit
this to set the target OSM infrastructure datasets, number of spatial slices, and
hazard datasets. See below and [config/README.md](https://github.com/nismod/open-gira/blob/main/config/README.md)
for more details on the configuration variables.

### Available pipelines

#### Network creation

`open-gira` can currently create several types of connected infrastructure
network from open data.

##### Road

We can create a topologically connected road network for a given area from
[OpenStreetMap](https://www.openstreetmap.org) (OSM) data. The resulting
network can be annotated with data retrieved from OSM (e.g. highway
classification, surface type), along with data looked up from user-supplied
sources (e.g. rehabilitation costs). The network edges will be labelled with
from nodes and to nodes, describing the connectedness of the network.

To specify a desired network:
- Review and amend the spreadsheets in `bundled_data/transport`, these supply
  information that is used to gap-fill or extend what can be determined from OSM alone.
- Review and amend `config/config.yaml`:
    - The `infrastructure_datasets` map should contain a key pointing to an `.osm.pbf`
      file URL for desired area. There are currently entries for the planet,
      for (some definition of) continents and several countries. We use
      the [geofabrik](http://download.geofabrik.de/) service for continent and
      country-level OSM extracts.
    - Check the OSM filter file pointed to by `network_filters.road`.
      This file specifies which [elements](https://wiki.openstreetmap.org/wiki/Elements)
      (nodes, ways or relations) to keep (or reject) from the multitude of data
      in an OSM file. See the filter expressions section
      [here](https://docs.osmcode.org/osmium/latest/osmium-tags-filter.html)
      for more information on the syntax of these files.
    - Check and amend `keep_tags.road`. This list of strings specifies which
      `tags` (attributes) to retain on the filtered elements we extract from
      the `.osm.pbf` file.
    - Review `slice_count`. This controls the degree of parallelism possible.
      With it set to 1, there is no spatial slicing (we create the network in
      a single chunk). To speed network creation for large domains, it can be
      set to a larger square number. The first square number greater than your
      number of available CPUs is a good heuristic.
    - Check and amend the values of `transport.road`, which provide some
      defaults for OSM data gap-filling.

And to create the network, by way of example:
```
snakemake --cores all -- results/egypt-latest_filter-road/edges.geoparquet
```

##### Rail

The process for creating a rail network is essentially the same as for road.
Please see the road section above for the relevant options that can be configured.

An example network creation call would be:
```
snakemake --cores all -- results/egypt-latest_filter-rail/edges.geoparquet
```

Note that the nodes file, `results/egypt-latest_filter-rail/nodes.geoparquet`
will by default contain the stations and their names.

##### Electricity grid creation

TODO

#### Risk assessment

TODO

##### Transport / flooding

The pipeline starts from a OpenStreetMap dataset (_e.g._
`europe-latest`) and produces network/flood hazard intersection data,
associating road splits to corresponding flood levels.

The pipeline consists in the following steps:

1. The target OSM datasets are downloaded or copied and saved as
   `<output_dir>/input/<dataset>.osm.pbf`.
2. The initial OSM datasets are filtered, keeping only relevant tags for road links
   (using `osmium tags-filter`). This results in smaller files
   `<output_dir>/input/<dataset>_filter-<filters>.osm.pbf`, where `<dataset>` is the
   key name and `<filters>` is the filename of the `osmium_tags_filter` file in the config.
3. The OSM dataset's headers are examined for a `bbox` property and that is used
   to determine the bounding box for the whole area (`<output_dir>/json/<dataset>.json`).
4. The hazard raster files for each hazard datasets are located by reading a list
   of their locations from the config. Each of these locations is visited and the
   .tif file downloaded or copied to `<output_dir>/input/hazard-<hazard>/raw/<filename>`
   where `<hazard>` is the keyname in the config and `<filename>` is the file's
   base name.
5. Each hazard raster file is clipped to contain just the hazard data for each dataset.
   These files are stored in `<output_dir>/input/hazard-<hazard>/<dataset>/<filename>`
   where `<dataset>` is the OSM dataset whose bounding box is used for clipping.
6. The OSM dataset bounding box is sliced into a grid of smaller bounding boxes
   according to the `slice_count` config option, and these slices are saved
   in a json file `<output_dir>/json/<dataset>-extracts.geojson`.
7. The filtered OSM file is sliced into areas of equal size using the bounding
   box grid from step 6. The slices are saved to
   `<output_dir>/slices/<dataset>_filter-<filter>/slice-<N>.osm.pbf`.
8. Each filtered OSM dataset slice is then converted to the GeoParquet data format,
   resulting in `<output_dir>/geoparquet/<dataset>_filter-<filters>_slice-<N>.geoparquet`.
9. Each geoparquet slice is intersected against flood level data from the
   hazard datasets. The hazard datasets consist of a collection of
   raster data files. The network/hazard intersection results in data
   `<output_dir>/splits/<dataset>_filter-<filters>_slice-<N>_hazard-<hazard>.geoparquet`
   describing roads split according to the raster grid and associated flood level values.
   A corresponding `parquet` files (without geometries) is also created.
10. Split data is then joined into a unique dataset describing
    infrastructure and associated hazard level values for each combination of
    OSM dataset and hazard dataset. This results in
    `<output_dir>/<dataset>_filter-<filters>_hazard-<hazard>.geoparquet`.
11. Each .geoparquet file is processed to produce a .tif raster file showing the length
    of road affected by flooding greater than a threshold defined in the config.
    These files are saved as `<output_dir>/exposure/<dataset>_filter-<filters>/hazard-<hazard>/raster/exposure_<hazard_tif_filename>`
12. Coastline data are downloaded from the NaturalEarthData.com service and saved in
    `<output_dir>/input/coastlines/ne_10m_ocean`.
13. Administrative boundary data are downloaded from the NaturalEarthData.com and saved in
    `<output_dir>/input/admin-boundaries/zip/` and extracted to `<output_dir>/input/admin-boundaries/ne_50m/`
14. The raster file from step 11 is combined with the coastline and admin boundary data from steps 12 and 13,
    to produce an overall image showing the raster data in its geographical context, located in
    `<output_dir>/exposure/<dataset>_filter-<filters>/hazard-<hazard>/img/exposure_<hazard_tif_base_filename>.png`

This is a directional acyclic graph (DAG) of a simplified version of the workflow
that uses just one OSM dataset, one hazard dataset, and one slice:
![DAG of the Snakefile workflow](docs/src/img/DAG-simple.png)

##### Electricity grid / tropical cyclone

TODO

#### Cleaning intermediate outputs

You can remove intermediate files by running the `clean` rule. To check what will be deleted,
```bash
snakemake -c1 -R clean
```

### Utilities

`open-gira` comes with a few small utilities outside the `snakemake` workflows.

#### Geoparquet -> Geopackage

As standard we use the `.geoparquet` format to store vector data on disk.
Unfortunately common GIS software such as QGIS may not yet support this file
format. To convert file(s) to geopackage, use:
```
python workflow/scripts/pq_to_gpkg.py <path_to_geoparquet_1> <path_to_geoparquet_2> <...>
```

This will write `.geopackage` files beside their source `.geoparquet`.

#### Unpickling interactive plots

`matplotlib` plots can be interactive (zoom, pan, etc.), but not as static
images. Some rules produce pickled plot files. To view these, use:
```
python workflow/scripts/unpickle_plot.py <path_to_pickled_plot>
```

## Documentation

Documentation is written using the [`mdbook`](https://rust-lang.github.io/mdBook/index.html)
format, using markdown files in the `./docs` directory.

Follow the [installation instructions](https://rust-lang.github.io/mdBook/guide/installation.html)
to get the `mdbook` command-line tool.

Build the docs locally:

```bash
cd docs
mdbook build
open book/index.html
```

Or run `mdbook serve` to run a server and rebuild the docs as you make changes.

## Related projects

Two libraries have been developed in tandem with `open-gira` and provide some
key functionality.

### snail

The open-source Python library [snail](https://github.com/nismod/snail)
is used for vector-raster intersection, e.g. identifying which road segments
might be affected by a set of flood map hazard rasters.

### snkit

The [snkit](https://github.com/tomalrussell/snkit) library is used for
network cleaning and assembly.

## Acknowledgments

This research received funding from the FCDO Climate Compatible Growth
Programme. The views expressed here do not necessarily reflect the UK
government's official policies.
