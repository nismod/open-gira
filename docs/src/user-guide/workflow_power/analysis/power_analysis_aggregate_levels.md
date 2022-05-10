# Statistic aggregation



For context, see the [target statistic gathering](power_analysis_target_specific.md) workflow first.
Based on the `results/power_output/statistics/aggregate/targets_geo_top{top_percentile_select}{T or F}percent.gpkg` file, this script now generates a new gpkg file: `results/power_output/statistics/aggregate/targets_geo_top{top_percentile_select}{T or F}percent_aggregated_region.gpkg`
which has level geometries with the aggregated statistics as features. The level must be specified in the config.yaml file
as the `aggregate_level` parameter. See details on the level selection in the [download levels file](../download/power_download_adminboundaries.md) and/or https://gadm.org/download_world.html. 
The parameter selection `aggregate_level = 0` equates to countries and so should a finder region be requested, then  `aggregate_level = 1` is recommended.
 `aggregate_level = 2` and above is feasible but comes at the cost of a longer computational time. An example for `aggregate_level = 1` is seen below (note that the legend must be generated through QGIS or elsewise).

![Aggregate statistics with `aggregate_level = 1` example](../power_img/aggregate.png)


The entire gathering and aggregation workflow can be called through the `analyse_aggregate_levels` rule.
