"""
The variables in this file are used throughout the power analysis rules.

If we can eventually do without them entirely, that would be great.
"""


from typing import List, Tuple


#### POWER/STORMS WORKFLOW ####

# how many samples is each storm track dataset split into?
SAMPLES_PER_TRACKSET = {
    "IBTrACS": 1,
    "STORM": 10,
    "IRIS": 10,
}
