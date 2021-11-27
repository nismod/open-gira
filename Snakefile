from glob import glob
# ------
# Read directories from config file
configfile: "config.yaml"


DATA_DIR = config["data_dir"]
OUTPUT_DIR = config["output_dir"]
HAZARD_DATA_DIR = config["hazard_data_dir"]
DATASET = config["dataset"]
hazard_slug = os.path.basename(config["hazard_csv"]).replace(".csv", "")

rule all:
    input:
        os.path.join(OUTPUT_DIR, f"{DATASET}.highway-core_{hazard_slug}_splits.geoparquet")

include: "rules/slice.smk"
include: "rules/filter_osm_data.smk"
include: "rules/convert_to_geoparquet.smk"
include: "rules/network_hazard_intersection.smk"
include: "rules/join_data.smk"

rule clean:
    shell:
        """
        rm -rf data/slices &&
        rm -rf outputs/slices
        """
