# Workflow

All of the work that `open-gira` does is encoded in a workflow that can be run by the tool
[`snakemake`](https://snakemake.readthedocs.io/en/stable/).

The workflow consists of four main steps:
1. [Download the data](download.md)
2. [Process the downloaded data](process.md)
3. [Intersect the hazards and infrastructure](intersect.md)
4. [Analyse the intersection](analysis.md)

There are currently two main strands of analysis, which may be reorganised in the future to
be more consistent and extended to include other infrastructure networks and other hazards.

## Transport-flooding analysis

The steps in the workflow process:
1. Download OpenStreetMap data
2. Filter OpenStreetMap data to focus on major infrastructure components
3. Determine the bounding box from the OpenStreetMap data
4. Download hazard raster files
5. Clip hazard raster files to bounding boxes determined in step 3
6. Calculate a grid of bounding boxes for slicing the OpenStreetMap data
7. Slice the OpenStreetMap data into smaller sections
8. Convert OpenStreetMap data to .geoparquet format
9. Add hazard information to infrastructure geoparquet
10. Join slices together to produce overall geoparquet file

These steps, along with the output produced at each stage,
are described in the subsections of this chapter.

These steps are summarised in the digital acyclic graph for `slice_count: 1`, for just the
`tanzania-latest` infrastructure and `aqueduct-coast` hazard data: [![DAG of the workflow for
the Tanzania dataset and coast flooding data](./img/DAG-simple.png)](./img/DAG-simple.png)

## Power-cyclones analysis

The power analysis workflow runs from the downloading of the data to the analysis of the
intersection. First the parameters need to be set in the config.yaml file
(`open-gira/config/config.yaml`). Each parameter is explained in this file and further
explanations are found in the relevant workflow documentation. All documentation files for the
power analysis will assume that `output_dir = results` when specifying paths.

The workflow for the power analysis consists of 4 main steps. Each of these steps can be called
through the linux command line when in the `open-gira` directory:

```shell
snakemake -s workflow/Snakefile {rule command} -{core numbers}
```

The `{rule command}` options are as follows
1. Download: `download_all`
2. Process: `process_all`
3. Intersect: `intersect_all`
4. Analyse: `analyse_all`

Note that the `figures_all` command to generate the figures (c.f. [here](analysis/figures.md))
can only be used once the main 4 steps (see above) have been completed.

By running any command, the previous ones will be included (e.g. if `intersect_all` is called,
then `download_all` and `process_all` are also performed). The `{core_numbers}` represents the
maximum number of cores that snakemake will be allocated. Furthermore, adding `-n` at the end
will allow for a 'dry-run' in which the workflow processes are shown but not executed. This
allows for the user to check expected workflow outputs. Lastly, should the workflow be
interruped (not recommended), then `--rerun-incomplete` should be added on the end of the
command.

If issues occur, it is recommended to do each of the four commands separately. Each workflow
process has a further detailed explanation and possible issue workarounds.
