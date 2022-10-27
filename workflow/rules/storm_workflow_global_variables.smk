"""
The variables in this file are used throughout the power analysis rules.

If we can eventually do without them entirely, that would be great.
"""


import requests
from typing import List, Tuple


def country_codes() -> List[str]:
    """
    Get a list of ISO A3 country codes from worldpop.org

    Returns:
        list[str]
    """

    # scripted requests are sometimes responded to with 403 permission denied
    # changing to one from a browser will circumvent the access control and return 200
    headers = {
        'user-agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:100.0) Gecko/20100101 Firefox/100.0',
    }
    r = requests.get("https://www.worldpop.org/rest/data/pop/cic2020_UNadj_100m", headers=headers)

    return [row["iso3"] for row in r.json()["data"]]


def all_boxes() -> List[str]:
    """
    Generate a list of box IDs for the cyclones workflow

    Returns:
        list[str]
    """
    return [f"box_{num}" for num in all_box_ids()]


def all_box_ids() -> List[int]:
    if len(config["specific_boxes"]) != 0:
        return config["specific_boxes"]

    max_i = int(360 * 180 / float(config["box_width_height"]) ** 2)
    return range(0, max_i)


#### POWER/STORMS WORKFLOW ####

# list of ISO A3 country codes
COUNTRY_CODES = country_codes()

# list of IDs of form "box_<int>"
ALL_BOXES = all_boxes()

CONNECTOR_OUT = (
    expand(
        os.path.join(
            config["output_dir"],
            "processed",
            "power",
            "{box_id}",
            "connector_{box_id}.json",
        ),
        box_id=all_box_ids(),
    )
)

STORM_BASINS = config["storm_basins"]
if len(STORM_BASINS) == 0:
    print("Inputting all storm basins")
    # east pacific, north atlantic, north indian, south india, south pacific, west pacific
    STORM_BASINS =  ("EP", "NA", "NI", "SI", "SP", "WP")

SAMPLES = config["storm_files_sample_set"]
if not SAMPLES:
    # empty list interpreted as 'run with all available samples'
    SAMPLES = list(range(0, 10))

STORMS = config["specific_storm_analysis"]
if STORMS == "None":
    STORMS = None

STORM_RPS = (
    list(range(10, 100, 10))
    + list(range(100, 1000, 100))
    + list(range(1000, 11000, 1000))
)

STORM_GCMS = ("CMCC-CM2-VHR4", "CNRM-CM6-1-HR", "EC-Earth3P-HR", "HadGEM3-GC31-HM")

STORM_RETURN_PERIODS_CURRENT = expand(
    os.path.join(
        config["output_dir"],
        "input",
        "storm-ibtracs",
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
        "storm-ibtracs",
        "fixed",
        "{storm_gcm}",
        "STORM_FIXED_RETURN_PERIODS_{storm_gcm}_{storm_rp}_YR_RP.tif",
    ),
    storm_gcm=STORM_GCMS,
    storm_rp=STORM_RPS,
)

STORM_RETURN_PERIODS = STORM_RETURN_PERIODS_CURRENT + STORM_RETURN_PERIODS_FUTURE

STORM_EVENTS_CURRENT = expand(
    os.path.join(
        config["output_dir"],
        "input",
        "storm-ibtracs",
        "events",
        "constant",
        "{storm_basin}",
        "STORM_DATA_IBTRACS_{storm_basin}_1000_YEARS_{num}.txt",
    ),
    storm_basin=STORM_BASINS,
    num=SAMPLES,
)

STORM_EVENTS_FUTURE = expand(
    os.path.join(
        config["output_dir"],
        "input",
        "storm-ibtracs",
        "events",
        "{storm_gcm}",
        "{storm_basin}",
        "STORM_DATA_{storm_gcm}_{storm_basin}_1000_YEARS_{num}_IBTRACSDELTA.txt",
    ),
    storm_gcm=STORM_GCMS,
    storm_basin=STORM_BASINS,
    num=SAMPLES,
)

STORM_EVENTS = STORM_EVENTS_CURRENT + STORM_EVENTS_FUTURE

try:
    STORM_BATCH_SIZE = int(config["storm_batches"])
    assert STORM_BATCH_SIZE > 0
except:
    raise RuntimeError("storm_batches incorrectly specified in config.yaml file")

WIND_RERUN_BOOL = config["wind_rerun"]
assert WIND_RERUN_BOOL in [True, False]

def get_storm_file(wildcards):
    """Helper to get storm events file, given:
    - OUTPUT_DIR
    - STORM_BASIN
    - STORM_MODEL (global climate model)
    - SAMPLE (0-9)
    """
    if wildcards.STORM_MODEL == "constant":
        fname = f"{wildcards.OUTPUT_DIR}/input/storm-ibtracs/events/constant/{wildcards.STORM_BASIN}/STORM_DATA_IBTRACS_{wildcards.STORM_BASIN}_1000_YEARS_{wildcards.SAMPLE}.txt"
    else:
        fname = f"{wildcards.OUTPUT_DIR}/input/storm-ibtracs/events/{wildcards.STORM_MODEL}/{wildcards.STORM_BASIN}/STORM_DATA_{wildcards.STORM_MODEL}_{wildcards.STORM_BASIN}_1000_YEARS_{wildcards.SAMPLE}_IBTRACSDELTA.txt"
    return fname

# check wind speed thresholds for damage are correctly ordered
assert config["central_threshold"] >= config["minimum_threshold"]
assert config["central_threshold"] <= config["maximum_threshold"]
WIND_SPEED_THRESHOLDS_MS = [
    config["central_threshold"],
    config["minimum_threshold"],
    config["maximum_threshold"],
]

# these files are written by the storm intersection script on finishing
# they are essentially a flag indicating successful completion
COMPLETION_FLAG_FILES = expand(
    os.path.join(
        config["output_dir"],
        "power_intersection",
        "storm_data",
        "all_winds",
        "{storm_basin}",
        "{sample}",
        "completed.txt",
    ),
    sample=SAMPLES,
    storm_basin=STORM_BASINS,
)

STORM_STATS_BY_THRESHOLD = expand(
    os.path.join(
        config["output_dir"],
        "power_output",
        "statistics",
        "combined_storm_statistics_{thrval}.csv",
    ),
    thrval=WIND_SPEED_THRESHOLDS_MS,
)

STORM_STATS_BY_REGION_SAMPLE_THRESHOLD = expand(
    os.path.join(
        config["output_dir"],
        "power_output",
        "statistics",
        "{storm_basin}",
        "{sample}",
        "combined_storm_statistics_{storm_basin}_{sample}_{thrval}.csv",
    ),
    storm_basin=STORM_BASINS,
    sample=SAMPLES,
    thrval=WIND_SPEED_THRESHOLDS_MS,
)

STORM_IMPACT_STATISTICS_DIR = os.path.join(config["output_dir"], "power_output", "statistics")

# variables to analyse for each storm
STORM_ANALYSIS_METRICS = [
    "GDP losses",
    "targets with no power (f=0)",
    "population affected",
    "population with no power (f=0)",
    "effective population affected",
    "reconstruction cost",
]

# variables to analyse for targets (electricity consumers)
TARGET_ANALYSIS_METRICS = [
    "population_without_power",
    "effective_population",
    "affected_population",
    "mw_loss_storm",
    "f_value",
    "gdp_damage",
]

# bastardised string version of boolean: 'T' or 'F'
SORT_BY_INCREASING_SEVERITY = 'T' if (config["increased_severity_sort"] == True) else 'F'
