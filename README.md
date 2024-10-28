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

The repository comes with a `environment.yml` file describing the `conda` and
`PyPI` packages required to run `open-gira`. The `open-gira` developers
recommend using either [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html#micromamba)
or [mamba](https://mamba.readthedocs.io/en/latest/index.html) to install and
manage these `conda` packages.

Having installed one of the suggested package managers, to create the
`open-gira` conda environment:

```bash
micromamba create -f environment.yml -y
```

And to activate the environment:

```bash
micromamba activate open-gira
```

### Utilities

Some rules use the `wget` utility to download files.

On Linux or MacOS, you may already have the `wget` utility available. If not,
it should be possible to install with your usual package manager (e.g. apt,
MacPorts, brew), or else using micromamba:

```bash
micromamba install wget
```

On Windows, you may have it already if you have a MinGW or Cygwin installation.
If not, you can access binaries at [eternallybored.org](https://eternallybored.org/misc/wget/).
Download the standalone exe and place it for example in `C:\Users\username\bin`
or somewhere on your PATH.

`exactextract` is used for zonal statistics in the tropical cyclones /
electricity grid analysis. It is not available via the `conda` package
management ecosystem and so must be installed separately. Please see
installation instructions [here](https://github.com/isciences/exactextract).

You are now ready to request result files, triggering analysis jobs in the
process.

Note that all subsequent commands given in the documentation assume that the
`open-gira` environment is already activated.

## Tests

Workflow steps are tested using small sample datasets.

To run the tests:

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
snakemake --cores 2 -- results/wales-latest_filter-road-primary/edges.gpq
```

Here, we ask `snakemake` to use up to 2 CPUs to produce a target file, in this
case, the edges of the Welsh road network. `snakemake` pattern matches
`wales-latest` as the OSM dataset name and `road-primary` as the network
type we want to filter for, picking up the [filter expressions](https://docs.osmcode.org/osmium/latest/osmium-tags-filter.html#filter-expressions) as defined in `config/osm_filters/road-primary.txt`.

To check what work we're going to request before commencing, use the `-n` flag:

```bash
snakemake -n --cores 2 -- results/wales-latest_filter-road-primary/edges.gpq
```

This will explain which rules will be required to run to produce the target
file. It may be helpful to [visualise](https://snakemake.readthedocs.io/en/stable/executing/cli.html#visualization)
which rules are expected to run, too.

The workflow configuration details are in `config/config.yml`. You can edit
this to set the target OSM infrastructure datasets, number of spatial slices, and
hazard datasets.

See the [documentation](https://nismod.github.io/open-gira/)
and [config/README.md](https://github.com/nismod/open-gira/blob/main/config/README.md)
for more details on usage in general and on configuration.

## Documentation

Documentation is written using the [`mdbook`](https://rust-lang.github.io/mdBook/index.html)
format, using markdown files in the `./docs` directory.

Follow the [installation instructions](https://rust-lang.github.io/mdBook/guide/installation.html)
to get the `mdbook` command-line tool.

To build the docs locally:

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

This research has also been supported by funding from the World Bank
Group, and the UK Natural Environment Research Council (NERC) through
the UK Centre for Greening Finance and Investment (CGFI).
