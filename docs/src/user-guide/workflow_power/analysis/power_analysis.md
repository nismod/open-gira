# Analysis of Results

Now that the downloading, processing and intersection workflow steps have been completed we can proceed with the analysis.
The following options are available

- [Empirical distribution generation](user-guide/workflow_power/analysis/power_analysis_empirical_distribution.md)
- [Empirical affected country pair matrix](user-guide/workflow_power/analysis/power_analysis_empirical_matrix.md)
- [Target statistic gathering](user-guide/workflow_power/analysis/power_analysis_target_specific.md) followed by [Statistic aggregation](user-guide/workflow_power/analysis/power_analysis_aggregate_levels.md)

These steps can be executed in one by entering `snakemake -s workflow/Snakefile analyse_all -c{core_numbers}`. Note
that if the `intersect_all`, `process_all` and/or `download_all` rules have not already run, then they will be run through this command
automatically.

