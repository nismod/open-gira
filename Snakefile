from glob import glob
# ------
# Read directories from config file
configfile: "config.yaml"


DATA_DIR = config["data_dir"]
OUTPUT_DIR = config["output_dir"]
HAZARD_DATA_DIR = config["hazard_data_dir"]
DATASET = config["dataset"]

# Variables for pattern rules

# For instance a full (unfiltered) osm dataset is named like {slug}
# where {slug} is e.g. tanzania-latest or
# tanzania-latest.highway-core. Defining the general structure of
# names using wildcards is useful to write general rules instead of
# hardcoding names. See
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards
hazard_slug = os.path.basename(config["hazard_csv"]).replace(".csv", "")

# Initial and final input file

INPUT_FILE = os.path.join(DATA_DIR, f"{DATASET}.osm.pbf")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, f"{DATASET}.highway-core_{hazard_slug}_splits.geoparquet")
INPUT_JSON_FILE = INPUT_FILE.replace(".osm.pbf", ".geojson")
EXTRACTS_CONFIG_FILE = INPUT_FILE.replace(".osm.pbf", "-extracts.geojson"),


rule all:
    input:
        OUTPUT_FILE,

include: "rules/slice.smk"
include: "rules/filter_osm_data.smk"
include: "rules/convert_to_geoparquet.smk"
include: "rules/network_hazard_intersection.smk"



def aggregate_input_geoparquet(wildcards):
    checkpoint_output = checkpoints.slice.get(**wildcards).output[0]
    return expand(
        os.path.join(
            OUTPUT_DIR,
            "slices",
            f"{DATASET}-slice{{i}}.highway-core_{hazard_slug}_splits.geoparquet",
        ),
        i=glob_wildcards(os.path.join(checkpoint_output, f"{DATASET}-slice{{i,\d+}}.osm.pbf")).i
    )

rule join_data:
    input:
        aggregate_input_geoparquet
    output:
        OUTPUT_FILE,
    script:
        "join_data.py"


rule clean:
    shell:
        """
        rm -rf data/slices &&
        rm -rf outputs/slices
        """
