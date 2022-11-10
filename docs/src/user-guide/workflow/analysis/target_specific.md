# Target statistic gathering


The target statistic gathering process iterates
through all targets and gathers (for each target) the sum and mean of each metric. In the config.yaml file, the
parameter `top_percentile_select` defines the percentile to analyse (in percent) for geospatial analysis and visualisation. Note
that `top_percentile_select = 100` will include all results. Also note that this percentile includes all storms even if they
have not met the threshold to be documented in the [wind analysis](../intersect/windspeed.md). It is recommended to
leave `increased_severity_sort = True`, however, should a top fraction of storms wish to be investigated and analysed only
then `increased_severity_sort = False` will reverse the percentile selection.

A new gpkg file is generated with the targets and their summed and average parameters as a target feature. This file is found
as `results/power_output/statistics/aggregate/targets_geo_top{top_percentile_select}{T or F}percent.gpkg` where `T` is selected
if `increased_severity_sort = True` and `F` if `increased_severity_sort = False` to ensure no data is overwritten.

One option now to visualise this data is to use the aggregate function in QGIS (or equivalent), however, this is a very time intensive process
which is not found to be very flexible. Rather, it is recommended to use the [statistic aggregation](aggregate_levels.md)
script. Click [here](aggregate_levels.md) to be directed onwards including
details on the rule to call the workflow.
