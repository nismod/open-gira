# Empirical distribution generation

The empirical distribution generation takes
the storm data and plots (according to the metric) the empirical distribution. The metrics which are analysed are
- GDP losses
- targets with no power (f=0)
- population affected
- population with no power (f=0)
- effective population affected (EACA see paper)
- Reconstruction costs (EAD see paper)

An example for EACA and EAD are seen (sequentially) below.
![EACA example](../../img/EACA.png)
![EAD example](../../img/EAD.png)


Each metric is saved to an image plot in `results/power_output/statistics/empirical`.

These files can be generated (in one go) through the `analyse_empirical_distribution` rule.
