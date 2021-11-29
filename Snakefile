from glob import glob

configfile: "config.yaml"

DATA_DIR = config["data_dir"]
OUTPUT_DIR = config["output_dir"]
HAZARD_DATA_DIR = config["hazard_data_dir"]
DATASET = config["dataset"]
hazard_slug = os.path.basename(config["hazard_csv"]).replace(".csv", "")

##### load rules #####

include: "workflow/rules/slice.smk"
include: "workflow/rules/filter_osm_data.smk"
include: "workflow/rules/convert_to_geoparquet.smk"
include: "workflow/rules/network_hazard_intersection.smk"
include: "workflow/rules/join_data.smk"

##### target rules #####

rule all:
    input:
        os.path.join(OUTPUT_DIR, f"{DATASET}.highway-core_{hazard_slug}_splits.geoparquet")

rule clean:
    shell:
        """
        rm -rf data/slices &&
        rm -rf outputs/slices
        """
