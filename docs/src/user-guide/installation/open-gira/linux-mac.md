# Linux/Mac installation

The major installation steps are to:
1. Download `open-gira`
1. Set up a Python environment
1. Install additional command-line tools (in particular, `osmium`)

## open-gira

Install open-gira by cloning the repository:

```bash
git clone https://github.com/nismod/open-gira.git
```

## Python environment

This repository comes with a `environment.yml` file describing the Python packages required to
run `open-gira`.

There are several ways to manage Python versions and install libraries.
- First tutorial and introduction to [installing Python
  packages](https://packaging.python.org/en/latest/tutorials/installing-packages/)
- [`conda`](https://docs.conda.io/en/latest/) lets you install different versions of Python
  and Python libraries and other dependencies.
- [`mamba`](https://mamba.readthedocs.io/en/latest/) is a replacement for `conda` which aims
  to do the same thing, faster.

The recommended approach for `open-gira` is to [install
`mamba`](https://mamba.readthedocs.io/en/latest/installation.html) then use it to create and
manage environments

Create the `open-gira` conda environment:

```bash
mamba env create -f environment.yml
```

and activate it

```bash
conda activate open-gira  # note that you still use `conda` to activate and deactivate
```

## Command-line tools

Install [`osmium`](https://osmcode.org/osmium-tool/manual.html) according to the
instructions there, or see the
[`Dockerfile`](https://github.com/nismod/open-gira/blob/main/Dockerfile) for example commands
to run on Ubuntu.

Install [`gdal` binaries](https://gdal.org/download.html)
- via conda/mamba, e.g. `conda install -c conda-forge gdal` should work cross-platform
- via system package manager, e.g. `sudo apt install gdal-bin` on Ubuntu.
