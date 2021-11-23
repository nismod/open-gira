# ------
# Read directories from config file
configfile: "config.yaml"


DATA_DIR = config["data_dir"]
OUTPUT_DIR = config["output_dir"]
AQUEDUCT_DIR = config["aqueduct_dir"]
DATASET = config["dataset"]

# ------
# Define ratio for slicing osm data into smaller areas
RATIO = config['ratio']
NSLICES = RATIO * RATIO

# ------
# Define naming of inputs/outputs at various stages of the
# pipeline
ALL_SLICE_FILES = [
    os.path.join(DATA_DIR, f"{DATASET}-slice{s}.osm.pbf")
    for s in range(0, NSLICES)
]
hazard_slug = os.path.basename(config["hazard_csv"]).replace(".csv", "")
ALL_GEOPARQUET_SPLITS_FILES = [
    slice_filename.replace(".osm.pbf", f".highway-core_{hazard_slug}_splits.geoparquet").replace(DATA_DIR, OUTPUT_DIR)
    for slice_filename in ALL_SLICE_FILES
]
ALL_PARQUET_SPLITS_FILES = [
    slice_filename.replace(".osm.pbf", f".highway-core_{hazard_slug}_splits.parquet").replace(DATA_DIR, OUTPUT_DIR)
    for slice_filename in ALL_SLICE_FILES
]

# Variables for pattern rules

# For instance a full (unfiltered) osm dataset is named like {slug}
# where {slug} is e.g. tanzania-latest or
# tanzania-latest.highway-core. Defining the general structure of
# names using wildcards is useful to write general rules instead of
# hardcoding names. See
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards
FULL_PBF_FILE = os.path.join(DATA_DIR, "{slug}.osm.pbf")
PBF_FILE = os.path.join(DATA_DIR, "{slug}.highway-core.osm.pbf")
GEOPARQUET_FILE = PBF_FILE.replace(".osm.pbf", ".geoparquet")
GEOPARQUET_SPLITS_FILE = GEOPARQUET_FILE.replace(
    ".geoparquet", f"_{hazard_slug}_splits.geoparquet"
).replace(DATA_DIR, OUTPUT_DIR)
PARQUET_SPLITS_FILE = GEOPARQUET_SPLITS_FILE.replace(".geoparquet", ".parquet")

# Initial and final input file

INPUT_FILE = os.path.join(DATA_DIR, f"{DATASET}.osm.pbf")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, f"{DATASET}.highway-core_{hazard_slug}_splits.geoparquet")
INPUT_JSON_FILE = INPUT_FILE.replace(".osm.pbf", ".geojson")
EXTRACTS_CONFIG_FILE = INPUT_FILE.replace(".osm.pbf", "-extracts.geojson"),


rule all:
    input:
        OUTPUT_FILE,

rule slice:
    input:
        data=INPUT_FILE,
        config=EXTRACTS_CONFIG_FILE
    output:
        ALL_SLICE_FILES,
    shell:
        "osmium extract --no-progress --config {input.config} {input.data}"


rule filter_osm_data:
    input:
        "filters.txt",
        FULL_PBF_FILE,
    output:
        PBF_FILE,
    shell:
        "osmium tags-filter {input[1]} w/highway=$(cat filters.txt) -o {output}"


rule convert_to_geoparquet:
    input:
        PBF_FILE,
    output:
        GEOPARQUET_FILE,
    script:
        "osm_to_pq.py"


rule network_hazard_intersection:
    input:
        network=GEOPARQUET_FILE,
    output:
        geoparquet=GEOPARQUET_SPLITS_FILE,
        parquet=PARQUET_SPLITS_FILE,
    script:
        "network_hazard_intersection.py"


rule join_data:
    input:
        data=ALL_GEOPARQUET_SPLITS_FILES,
    output:
        OUTPUT_FILE,
    script:
        "join_data.py"


rule clean:
    shell:
        """
        rm -f data/*-slice?.osm.pbf &&
        rm -f data/*-slice?.highway-core.osm.pbf &&
        rm -f data/*-slice?.highway-core.geoparquet &&
        rm -f outputs/*-slice?.highway-core.splits.geoparquet
        rm -f outputs/*-slice?.highway-core.splits.parquet
        """
