# Workflow for Power Analysis

The power analysis workflow runs from the downloading of the data to the analysis of the intersection. First the
parameters need to be set in the config.yaml file (`open-gira/config/config.yaml`). Each parameter is explained
in this file and further explanations are found in the relevant workflow documentation. All documentation files for the power
analysis will assume that `output_dir = results` when specifying paths.


The workflow for the power analysis consists of 4 main steps
1. [Download the data](download/power_download.md)
2. [Process the downloaded data](process/power_process.md)
3. [Intersect the downloaded storms with the processed power infrastructure](intersect/power_intersect.md)
4. [Analyse the intersection](analysis/power_analysis.md)

Each of these steps is explained in further detail in their respective links. They can be called through
the linux command line when in the `open-gira` directory:
```shell
snakemake -s workflow/Snakefile {rule command} -{core numbers}
```
The `{rule command}` options are as follows
1. Download: `download_all`
2. Process: `process_all`
3. Intersect: `intersect_all`
4. Analyse: `analyse_all`

By running any command, the previous ones will be included (e.g. if `intersect_all` is called, then `download_all` and `process_all` are also performed).
The `{core_numbers}` represents the maximum number of cores that snakemake will be allocated.
Furthermore, adding `-n` at the end will allow for a 'dry-run' in which the workflow processes are shown but not
executed. This allows for the user to check expected workflow outputs. Lastly, should the workflow be interruped
(not recommended), then `--rerun-incomplete` should be added on the end of the command.

If issues occur, it is recommended to doo each of the four commands separately. Each workflow process has a further detailed
explanation and possible issue workarounds. An overview of these files can be found [here](../../SUMMARY.md).

