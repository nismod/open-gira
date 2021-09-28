configfile: "config.yaml"

DATA_DIR = config["data_dir"]
OUTPUT_DIR = config["output_dir"]
AQUEDUCT_DIR = config["aqueduct_dir"]

FULL_PBF_FILE = os.path.join(DATA_DIR, "{slug}.osm.pbf")
PBF_FILE = os.path.join(DATA_DIR, "{slug}-highway-core.osm.pbf")
GEOPARQUET_FILE = PBF_FILE.replace(".osm.pbf", ".geoparquet")

GEOPARQUET_SPLITS_FILE = GEOPARQUET_FILE.replace(
    ".geoparquet", "_splits.geoparquet"
).replace(DATA_DIR, OUTPUT_DIR)

PARQUET_SPLITS_FILE = GEOPARQUET_SPLITS_FILE.replace(".geoparquet", ".parquet")


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


rule clean:
    shell:
        "rm -rf data/*-highway-core.osm.pbf data/*.geoparquet outputs/"
