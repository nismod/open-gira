# Chop hazard data according to the overall bounding box
# calculated in the .osm.pbf dataset's json file

from json import load

checkpoint trim_hazard_data:
    input:
        file="{OUTPUT_DIR}/input/{HAZARD_SLUG}/raw/{FILENAME}.tif",
        json=os.path.join(f"{config['output_dir']}", "json", f"{dataset_slug}.json")
    output:
        "{OUTPUT_DIR}/input/{HAZARD_SLUG}/{DATASET}/{FILENAME}.tif"
    run:
        out_file = os.path.join(
            f"{wildcards.OUTPUT_DIR}",
            "input",
            f"{wildcards.HAZARD_SLUG}",
            f"{dataset_slug}",
            f"{wildcards.FILENAME}.tif"
        )
        with open(input.json, "r") as j:
            dict = load(j)
            xmin, ymin, xmax, ymax = dict['extracts'][0]['bbox']
            os.system(f"gdalwarp -te {xmin} {ymin} {xmax} {ymax} {input.file} {out_file}")

# Not sure how to test this within snakemake.
# It can be done from the command line by requiring a specific file, e.g.:
# snakemake --cores all data/hazard-aqueduct-coast/tanzania-latest/inuncoast_historical_nosub_hist_rp0001_5.tif
