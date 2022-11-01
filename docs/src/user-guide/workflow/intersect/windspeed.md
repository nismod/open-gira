# Storm analysis

This file takes the [storm data](../download/storm-ibtracs.md) and calculates the wind speed parameters for
each unit (see [unit generation](grid.md)) for each storm.

The parameters investigated are
 - `wind_location` Wind speed at the unit at the storm time step
 - `wind` Wind at cyclone track at the storm time step
 - `wind_max` Maximum wind speed over entire storm at the cyclone track
 - `duration_15ms` Number of occurrences of wind speed above 15m/s at the unit
 - `duration_20ms` Number of occurrences of wind speed above 20m/s at the unit

Each time step is 3 hours. The following steps are applied:
1. load the [storm data](../download/storm-ibtracs.md)
2. For each [unit](grid.md) save the affected storms and their wind parameters to a parquet file
3. For each storm load the relevant units (see above), calculate the wind parameters at all (relevant) units and save to a .csv file

After this, each storm will have a .csv file containing all wind speed parameters at all relevant units. In addition, in step 2., only
storms with a storm centre is within `hurr_buffer_dist = 1300 km` (of the unit) will be saved. This can be visualised as only saving
data within a 'buffer' of the storm track. A unit that is too far from the storm centre will be unaffected by the storm. The
`hurr_buffer_dist` parameter is designed to capture this. Furthermore, in step 3., only storm parameters which have a
`wind_location` value above `min_windlocmax` and `wind_max` value above `min_windmax` are considered.
Values below this can be considered negligible and will not cause infrastructure damage.

Each storm has a unique identifier in the form #1_#2_#3 where #1 is the sample number from the [storm data](../download/storm-ibtracs.md), #2 is
the year in that sample and #3 is the #3-th storm of that year. This is often referred to by `nh` (number hurricane) within the workflow.


The processing time, memory and storage costs were heavily considered during the script development.As such measures have
been implemented to reduce these parameters. Rather than save the intermediate data to memory, each unit is saved to a
parquet file (step 2.) which can be later retrieved (step 3.). Furthermore, if unit parquet files already exist, it will
be skipped in step 2. and a dictionary is formed of which units contain which storm data in attempts to improve the
computational time.


Snakemake checkpoints (https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html) are used here as the number of
storms and their unique identifier is not known here. It is highly recommended not to interrupt this process as, in certain
cases, the checkpoint function will not know whether it has evaluated all storms and may not rerun in a following snakemake
command. Should this process be interrupted and you wish to rerun all wind evaluations, you must delete
`results/power_intersection/storm_data/all_winds/{region}/{sample}/{region}_{sample}_completed.txt` for the corresponding region(s) and sample(s).
