# MacOS installation

The most efficient way to install open-gira on MacOS is to use [conda](https://docs.conda.io/en/latest/).
Conda will take care of installing the tools on which open-gira relies.
If you do not wish to install via conda, you can install via pip instead and install the
osmium-tool libraries separately.

## conda

This repository comes with a `environment.yml` file describing the conda and pip packages required to run `open-gira`.

Create the `open-gira` conda environment:
```
conda env create -f workflow/envs/environment.yml
```
and activate it
```
conda activate open-gira
```

## pip

Pip installation requires installing the python requirements from `requirements.txt`,
as well as the osmium-tool software.

### Installing Python libraries

Install python requirements as listed in `requirements.txt` - for
example using a venv:

```
python3 -m venv ./venv
. venv/bin/activate
pip install -r requirements.txt
```

### Installing osmium-tool

Install [`osmium-tool`](https://osmcode.org/osmium-tool/manual.html) according to the
instructions there. Tests run with versions:
- osmium-tool v1.13.2
- libosmium v2.17.1

