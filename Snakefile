configfile: "config.yaml"

DATA_DIR = config["data_dir"]
OUTPUT_DIR = config["output_dir"]
AQUEDUCT_DIR = config["aqueduct_dir"]

FULL_PBF_FILE = os.path.join(DATA_DIR, "{slug}.osm.pbf")
PBF_FILE = os.path.join(DATA_DIR, "{slug}.highway-core.osm.pbf")
GEOPARQUET_FILE = PBF_FILE.replace(".osm.pbf", ".geoparquet")

GEOPARQUET_SPLITS_FILE = GEOPARQUET_FILE.replace(
    ".geoparquet", ".splits.geoparquet"
).replace(DATA_DIR, OUTPUT_DIR)

PARQUET_SPLITS_FILE = GEOPARQUET_SPLITS_FILE.replace(".geoparquet", ".parquet")

RATIO = 3
NSLICES = RATIO * RATIO

ALL_SPLITS_FILES = [
    os.path.join(
        OUTPUT_DIR, f"tanzania-latest-slice{s}.highway-core.splits.geoparquet"
    )
    for s in range(1, NSLICES+1)
]



rule all:
    input:
        os.path.join(OUTPUT_DIR, "tanzania-latest.splits.geoparquet")


rule slice:
    input:
        data="data/tanzania-latest.osm.pbf",
        cmd="split_to_bounding_boxes.sh"
    output: [f"data/tanzania-latest-slice{i}.osm.pbf" for i in range(1, NSLICES+1)]
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
    input: ALL_SPLITS_FILES
    output: os.path.join(OUTPUT_DIR, "tanzania-latest.splits.geoparquet")


rule clean:
    shell:
        """
        rm -f data/*-slice?.osm.pbf &&
        rm -f data/*-slice?.highway-core.osm.pbf &&
        rm -f data/*-slice?.highway-core.geoparquet &&
        rm -f outputs/*-slice?.highway-core.splits.geoparquet
        rm -f outputs/*-slice?.highway-core.splits.parquet
        """
