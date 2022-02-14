# Chop hazard data according to the overall bounding box
# calculated in the .osm.pbf dataset's json file

from json import load

# https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#data-dependent-conditional-execution
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.download_hazard_datasets.get(**wildcards).output[0]
    tifs = glob_wildcards(os.path.join(checkpoint_output, "{tif}.tif"))
    input = expand(os.path.join(checkpoint_output, "{tif}.tif"), tif=tifs.tif)
    # print(f"checkpoint_output={checkpoint_output}")
    # print(f"tifs={tifs}")
    # print(f"input={input}")
    return input

# This is a checkpoint because it outputs to a directory we want to ensure is up to date in later rules.
# Specifically, intersection.smk looks for all *.tif files in the output directory.
checkpoint trim_hazard_data:
    input:
        files=aggregate_input,
        json="{OUTPUT_DIR}/json/{DATASET}.json"
    output:
        directory("{OUTPUT_DIR}/input/{HAZARD_SLUG}/{DATASET}")
    run:
        os.system(f"mkdir --parents {output}")

        with open(input.json, "r") as j:
            dict = load(j)
            xmin, ymin, xmax, ymax = dict['extracts'][0]['bbox']

        for f in input.files:
            print(f"Trimming {f} for {wildcards.HAZARD_SLUG}, {wildcards.DATASET}.")
            out_file = f"{output}/{os.path.basename(f)}"
            os.system(f"gdalwarp -te {xmin} {ymin} {xmax} {ymax} {f} {out_file}")

"""
Test with:
snakemake --cores all results/input/hazard-aqueduct-river/tanzania-latest/inunriver_rcp4p5_MIROC-ESM-CHEM_2030_rp00100.tif
"""
