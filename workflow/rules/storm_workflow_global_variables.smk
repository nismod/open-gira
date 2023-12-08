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

STORM_RPS = (
    list(range(10, 100, 10))
    + list(range(100, 1000, 100))
    + list(range(1000, 11000, 1000))
)

STORM_GCMS = ("CMCC-CM2-VHR4", "CNRM-CM6-1-HR", "EC-Earth3P-HR", "HadGEM3-GC31-HM")
STORM_MODELS = STORM_GCMS + ("constant", )

STORM_RETURN_PERIODS_CURRENT = expand(
    os.path.join(
        config["output_dir"],
        "input",
        "STORM",
        "fixed",
        "constant",
        "STORM_FIXED_RETURN_PERIODS_constant_{storm_rp}_YR_RP.tif",
    ),
    storm_rp=STORM_RPS,
)

STORM_RETURN_PERIODS_FUTURE = expand(
    os.path.join(
        config["output_dir"],
        "input",
        "STORM",
        "fixed",
        "{storm_gcm}",
        "STORM_FIXED_RETURN_PERIODS_{storm_gcm}_{storm_rp}_YR_RP.tif",
    ),
    storm_gcm=STORM_GCMS,
    storm_rp=STORM_RPS,
)

STORM_RETURN_PERIODS = STORM_RETURN_PERIODS_CURRENT + STORM_RETURN_PERIODS_FUTURE
