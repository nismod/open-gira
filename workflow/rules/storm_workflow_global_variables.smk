"""
The variables in this file are used throughout the power analysis rules.

If we can eventually do without them entirely, that would be great, but this is
better than what came before.
"""


import requests
from typing import List


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

    if len(config["specific_boxes"]) != 0:
        return [f"box_{num}" for num in config["specific_boxes"]]

    else:
        return [
            f"box_{int(idx)}"
            for idx in range(
                0, int((180 - -180) * (90 - -90) / float(config["box_width_height"]) ** 2)
            )
        ]


#### POWER/STORMS WORKFLOW ####

# list of ISO A3 country codes
COUNTRY_CODES = country_codes()

# list of IDs of form "box_<int>"
ALL_BOXES = all_boxes()

ADMIN_BOUNDS_GLOBAL_SINGLE_LAYER = os.path.join(config['output_dir'], "input", "admin-boundaries", "gadm36.gpkg")
ADMIN_BOUNDS_GLOBAL_LAYER_PER_LEVEL = os.path.join(
    config['output_dir'], "input", "admin-boundaries", "gadm36_levels.gpkg"
)
ADMIN_BOUNDS_FILE_PER_COUNTRY = expand(
    os.path.join(
        config["output_dir"], "input", "admin-boundaries", "gadm36_{code}.gpkg"
    ),
    code=COUNTRY_CODES,
)

CONNECTOR_OUT = (
    expand(
        os.path.join(
            config["output_dir"],
            "power_processed",
            "all_boxes",
            "{box_id}",
            "connector_{box_id}.txt",
        ),
        box_id=ALL_BOXES,
    ),
)

# east pacific, north atlantic, north indian, south india, south pacific, west pacific
STORM_BASINS = ("EP", "NA", "NI", "SI", "SP", "WP")
REGIONS = config["regions"]
if len(REGIONS) == 0:
    print("Inputting all regions")
    REGIONS = STORM_BASINS

SAMPLES = list(range(config["sample_upper"] + 1))
if config["samples_indiv"] != "None":
    print("Using specified samples")
    SAMPLES = config["samples_indiv"]
if len(SAMPLES) == 0:
    raise ValueError("Samples incorrectly specified")

STORMS = config["specific_storm_analysis"]
if STORMS == "None":
    STORMS = None

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
        "{region}",
        "{sample}",
        "{region}_{sample}_completed.txt",
    ),
    sample=SAMPLES,
    region=REGIONS,
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
        "{region}",
        "{sample}",
        "combined_storm_statistics_{region}_{sample}_{thrval}.csv",
    ),
    region=REGIONS,
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
