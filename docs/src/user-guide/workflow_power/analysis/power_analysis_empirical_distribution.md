# Empirical distribution generation

The empirical distribution generation takes
the storm data and plots (according to the metric) the empirical distribution. The metrics which are analysed are 
- GDP losses
- targets with no power (f=0)
- population affected
- population with no power (f=0)
- effective population affected

An example for GDP losses and population affected are seen (sequentially) below.
![GDP Losses example](../power_img/empirical_GDP%20losses.png)
![Population affected example](../power_img/empirical_population%20affected.png)


Each metric is saved to an image plot in `results/power_output/statistics/empirical`.

These files can be generated (in one go) through the `analyse_empirical_distribution` rule.
