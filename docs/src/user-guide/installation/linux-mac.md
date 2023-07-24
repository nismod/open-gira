# Linux / Mac

The major installation steps are to:
1. Download `open-gira`
1. Set up a Python environment
1. Install additional command-line tools

## Clone repository

Install open-gira by cloning the repository:

```bash
git clone https://github.com/nismod/open-gira.git
```

## Software environment

This repository comes with a `environment.yml` file describing almost all of the
software dependencies required to run `open-gira`.

There are several ways to manage Python versions and install libraries.
- [`conda`](https://docs.conda.io/en/latest/) lets you install different versions of Python
  and Python libraries and other dependencies.
- [`micrombamba`](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html#micromamba) 
  a replacement for `conda`, and what the `open-gira` developers use.

The recommended approach for `open-gira` is to install `micromamba` then use it
to create and manage environments.

### Local

To install the required dependencies on a local machine, create the `open-gira`
conda environment:

```bash
micromamba create -f environment.yml -y
```

Then activate it:

```bash
micromamba activate open-gira
```

You're now ready to configure workflows and request outputs.

### Cluster

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

## Other command-line tools

The following tools are not available through `conda` and must be installed separately.

### exactextract

exactextract is used for zonal statistics. Please see installation instructions [here](https://github.com/isciences/exactextract).
