# ------
# Read directories from config file
configfile: "config.yaml"

DATA_DIR = config["data_dir"]
OUTPUT_DIR = config["output_dir"]
AQUEDUCT_DIR = config["aqueduct_dir"]
DATASET = config["dataset"]

# ------
# Define ratio for slicing osm data into smaller areas
RATIO = 3
NSLICES = RATIO * RATIO

# ------
# Define naming of inputs/outputs at various stages of the
# pipeline
ALL_SLICE_FILES = [
    os.path.join(
        OUTPUT_DIR, f"{DATASET}-slice{s}.osm.pbf"
    )
    for s in range(1, NSLICES+1)
]
ALL_GEOPARQUET_SPLITS_FILES = [slice_filename.replace(".osm.pbf", ".highway-core.geoparquet") for slice_filename in ALL_SLICE_FILES]
ALL_PARQUET_SPLITS_FILES = [slice_filename.replace(".osm.pbf", ".highway-core.parquet") for slice_filename in ALL_SLICE_FILES]

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
    ".geoparquet", ".splits.geoparquet"
).replace(DATA_DIR, OUTPUT_DIR)
PARQUET_SPLITS_FILE = GEOPARQUET_SPLITS_FILE.replace(".geoparquet", ".parquet")

# Initial and final input file

INPUT_FILE=os.path.join(DATA_DIR, f"{DATASET}.osm.pbf")
OUTPUT_FILE=os.path.join(OUTPUT_DIR, f"{DATASET}.splits.geoparquet")

rule all:
    input:
        OUTPUT_FILE


rule slice:
    input:
        data=INPUT_FILE,
        cmd="split_to_bounding_boxes.sh"
    output: ALL_SLICE_FILES
    shell: "bash {input.cmd} {input.data} {RATIO}"


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
        cmd="osm_to_pq.py",
        data=PBF_FILE,
    output:
        GEOPARQUET_FILE,
    shell:
        "python {input.cmd} {input.data} {DATA_DIR}"


rule network_hazard_intersection:
    input:
        cmd="network_hazard_intersection.py",
        network=GEOPARQUET_FILE,
        csv=os.path.join(AQUEDUCT_DIR, "aqueduct_river.csv"),
    output:
        geoparquet=GEOPARQUET_SPLITS_FILE,
        parquet=PARQUET_SPLITS_FILE,
    shell:
        "python {input.cmd} {input.network} {AQUEDUCT_DIR} {OUTPUT_DIR}"


rule join_data:
    input:
        data=ALL_GEOPARQUET_SPLITS_FILES,
        cmd="join_data.py"
    output: OUTPUT_FILE
    shell: "python {input.cmd} {input.data}"


rule clean:
    shell:
        """
        rm -f data/*-slice?.osm.pbf &&
        rm -f data/*-slice?.highway-core.osm.pbf &&
        rm -f data/*-slice?.highway-core.geoparquet &&
        rm -f outputs/*-slice?.highway-core.splits.geoparquet
        rm -f outputs/*-slice?.highway-core.splits.parquet
        """
