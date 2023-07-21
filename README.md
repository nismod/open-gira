# Open Global Infrastructure Risk/Resilience Analysis

[![mdBook Documentation](https://github.com/nismod/open-gira/actions/workflows/docs.yml/badge.svg?branch=main)](https://nismod.github.io/open-gira)
[![pyTest](https://github.com/nismod/open-gira/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/nismod/open-gira/actions/workflows/test.yml)
[![snakemake workflow](https://img.shields.io/badge/snakemake-open--gira-informational)](https://snakemake.github.io/snakemake-workflow-catalog/?usage=nismod/open-gira)

## Introduction

This open-source [snakemake](https://snakemake.readthedocs.io/en/stable/)
workflow can be used to analyse environmental risks to infrastructure
networks using global open data. It is a work in progress.

Goals:
- Automated pipeline for reproducible analysis anywhere in the world
- Maps per-country and of larger areas
- Charts/stats of exposure per admin region, per hazard type, scenario, epoch
- Consider transport, electricity, water, communications systems
- Consider river flooding, storm surge coastal flooding, tropical cyclones
- Estimate direct damages to physical networks
- Estimate indirect effects of disruption - people affected, economic activity disrupted

Non-goals:
- Using closed data, which may be appropriate for other projects or use-cases
- Detailed operational/engineering level simulation
- Long-term planning

## Installation

Install `open-gira` by cloning the repository:
```bash
git clone https://github.com/nismod/open-gira.git
```

### Conda packages

The repository comes with a `environment.yml` file describing the `conda` and
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

Note that all subsequent commands given in the documentation assume that the
`open-gira` environment is already activated.

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
```bash
snakemake --cores 2 -- results/wales-latest_filter-road/edges.geoparquet
```

Here, we ask `snakemake` to use up to 2 CPUs to produce a target file, in this
case, the edges of the Welsh road network. `snakemake` pattern matches
`wales-latest` as the OSM dataset name and `filter-road` as the network type we
want to filter for.

To check what work we're going to request before commencing, use the `-n` flag:
```bash
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

### Quick start

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
```bash
snakemake --cores all -- results/egypt-latest_filter-road/edges.geoparquet
```

##### Rail

The process for creating a rail network is essentially the same as for road.
Please see the road section above for the relevant options that can be configured.

An example network creation call would be:
```bash
snakemake --cores all -- results/egypt-latest_filter-rail/edges.geoparquet
```

Note that the nodes file, `results/egypt-latest_filter-rail/nodes.geoparquet`
will by default contain the stations and their names as recorded in OSM.

##### Electricity grid creation

To create an electricity grid we rely heavily on [gridfinder](https://gridfinder.rdrn.me/)
data. This dataset provides transmission and distribution edges with
substantial global coverage. It also contains a set of 'targets' or electricity
consuming areas, derived from [Night Time Lights](https://www.earthdata.nasa.gov/learn/backgrounders/nighttime-lights)
(NTL) satellite imagery. Our other major data source for electricity grid creation is the World Resources Institute's (WRI)
[powerplants database](https://datasets.wri.org/dataset/globalpowerplantdatabase).

The workflow currently divides network creation by country. One may request one
country, or several. Note that neighbouring countries' networks are _not_
connected with one another.

There aren't currently any options to configure when creating an electricity
grid.

Here's an example grid creation command for Puerto Rico:
```bash
snakemake --cores all -- results/power/by_country/PRI/network/edges.geoparquet
```

The folder name under `by_county` is an [ISO 3166 Alpha-3 country code](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3#Officially_assigned_code_elements),
specifying the country in question.

#### Risk assessment

`open-gira` is designed to create spatial representations of infrastructure
networks _and_ to analyse their disruption due to environmental hazards.
Currently the workflow supports studies of flooding (river, coastal) to road
and rail networks and also the effect of tropical cyclones on electricity
grids.

##### Transport / flooding

The flooding risk analysis pipeline starts by creating an infrastructure
network (road or rail) as described above. Please refer to this section to
configure the network creation.

The hazard component of the analysis is configurable in `config/config.yaml`:
- `hazard_datasets` contains hazard names pointing to files of hazard layers.
  These layers are currently flood rasters (inundation depths for a given return
  period). Add or amend an entry pointing to file containing the rasters you
  wish to consider.
- Ensure `hazard_types` contains an identical key referencing the hazard types.
  This is currently limited to `flood` only.
- Configure the damage curves:
    - Check and amend `direct_damages.asset_types` contains any assets you wish
      to calcuate direct damage costs for. Currently implemented assets are
      available in `src/open_gira/assets.py` as the classes inheriting from
      `Assets`.
    - Ensure `direct_damages.curves_dir` is set to a path containing damages
      functions per asset type, organised by hazard type. See
      `bundled_data/damage_curves` for an example.
    - These damage function files should be named `<asset_type>.csv`, e.g.
      `road_unpaved.csv`. Their format is exemplified by
      `bundled_data/damage_curves/flood/road_unpaved.csv`

To request an evaluation of Expected Annual Damages (EAD) as a function of
hazard Return Period (RP) for a given `slice`, we can request something like:
```bash
snakemake --cores all -- results/direct_damages/<dataset_name>_filter-<network_type>/hazard-<hazard_name>/EAD_and_cost_per_RP/slice-<slice_number>.geoparquet
```

For example (with a `config.slice_count` of 9):
```bash
snakemake --cores all -- results/direct_damages/egypt-latest_filter-road/hazard-aqueduct-river/EAD_and_cost_per_RP/slice-5.geoparquet
```

And to compute all the slices for a given domain and then aggregate to country level (admin level 0):
```bash
snakemake --cores all -- results/egypt-latest_filter-road/hazard-aqueduct-river/EAD_and_cost_per_RP/agg-sum/admin-level-0.geoparquet
```

For more possible outputs please refer to the detailed documentation and the
rules defined in `workflow/rules/`.

##### Electricity grid / tropical cyclone

This analysis intersects electricity grids with tropical cyclone wind speed
risk. Network creation is described above. The hazards are event-based, as
opposed to return period maps. We take a storm track, estimate the resulting
wind field as a function of time, and then take the maximum across time. This
produces a maximum wind speed field for a given event. For a certain defined
threshold or set of thresholds, we can then fail grid edges that experience
wind in excess of the threshold.

To configure the analysis:
- Review and amend `config/config.yaml`:
    - `storm_sets` should contain a storm set name, pointing to a JSON file.
      This JSON file should contain an empty list to process all storms for
      this storm set, or a list of string storm IDs if only a subset is
      required.
    - Specify `transmission_windspeed_failure` as a list of wind speeds in
      meters per second at which to fail edges.

See comments in the `config/config.yaml` file for other less crucial
configuration options.

To request an analysis for a given storm-country pair:
```bash
snakemake --cores 1 results/power/by_country/PRI/exposure/IBTrACS/2017242N16333.nc
```

Or for all countries affected by a storm:
```bash
snakemake -c1 results/power/by_storm_set/IBTrACS/by_storm/2017242N16333/exposure_by_target.nc
```

For more possible outputs please see the more detailed documentation.

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
