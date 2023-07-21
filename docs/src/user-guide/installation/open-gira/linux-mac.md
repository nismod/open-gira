# Linux/Mac installation

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
- First tutorial and introduction to [installing Python
  packages](https://packaging.python.org/en/latest/tutorials/installing-packages/)
- [`conda`](https://docs.conda.io/en/latest/) lets you install different versions of Python
  and Python libraries and other dependencies.
- [`mamba`](https://mamba.readthedocs.io/en/latest/) is a replacement for `conda` which aims
  to do the same thing, faster.
- [`micrombamba`](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html#micromamba) 
  another replacement for `conda`, and what the `open-gira` developers use.

The recommended approach for `open-gira` is to install `micromamba` then use it
to create and manage environments.

Create the `open-gira` conda environment:

```bash
micromamba create -f environment.yml -y
```

and activate it

```bash
micromamba activate open-gira
```

## Other command-line tools

The following tools are not available through `conda` and must be installed separately.

### exactextract

exactextract is used for zonal statistics. Please see installation instructions [here](https://github.com/isciences/exactextract).

