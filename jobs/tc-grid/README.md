# Tropical Cyclone electricity grid impact workflow

## Introduction

This set of scripts can be used to run an open-gira Tropical Cyclone (TC)
electricity grid impact analysis. This workflow creates models of electricity
grids in the global TC belt and many years worth of synthetic TC wind speed
footprints. The disruption of electricity supply due to these events is
estimated. Output quantities include the size of outage events (in people
affected) for given return periods for each country, along with the expected
annual people affected for a set of administrative areas.

## SLURM workflow

It is possible to orchestrate the whole workflow using snakemake alone,
creating the Directed Acyclic Graph (DAG) and running the required jobs subject
to their dependencies. This is best achieved on a large workstation (e.g. 64
CPU 128GB RAM will run 10,000 years of global tracks in ~1 day).

It is however difficult to run the workflow in this way on a SLURM cluster.
This is because the workflow uses dynamic resource allocation (some jobs need
to vary memory or CPUs based on input data), which does not currently (early
2026) work with the snakemake-SLURM plugin. As such, to run on SLURM HPC, we
split the workflow into three phases as a workaround.

### Preprocess
We first create the electricity grids, TC track subsets and wind downscaling
factors for all desired countries and storm sets.

This single SLURM job calls snakemake on a single node and runs the DAG to the
above point:
`sbatch jobs/tc-grid/0-preprocess.sh`

### Core
This phase calculates wind footprints for all TC events and then estimates the
disruption to the electricity network that results.

We use SLURM array jobs, one task for each country/millenium-sample/storm-set
combination. Each array task will copy the workflow rules and source code to a
new directory and run its own snakemake parent process to create the minimal
DAG required to produce the specified output file. The output files are written
to the one and only output directory (symlinked from each task's working
directory).

To generate the SLURM array submission scripts, run (see script for optional arguments):
`python jobs/tc-grid/1-core-generate-job-scripts.py`

This will generate submission scripts, split into three size categories (small,
medium and large). Each category has different resource requirements.

These can then be submitted to the SLURM scheduler. They do not have
interdependencies, so may be submitted simultaneously. However it is
recommended you test one of them with a small number of tasks first -- there
will potentially be many thousands of jobs!
`sbatch jobs/tc-grid/1-core-small.sh`
`sbatch jobs/tc-grid/1-core-medium.sh`
`sbatch jobs/tc-grid/1-core-large.sh`

### Postprocess

The last major phase of the analysis is to collate results within storm-sets
and compute return period and EAPA figures.

We use a single SLURM job to call snakemake on a single node and run the DAG
using the outputs from phase 1:
`sbatch jobs/tc-grid/2-postprocess.sh`

## Further analysis

The `open-gira` portion of the workflow is now complete. Subsequent analysis of
the results (plotting of figures for papers etc.) is stored elesewhere.
