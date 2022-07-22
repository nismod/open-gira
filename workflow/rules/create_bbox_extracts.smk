# Create *-extracts.json file that determines the bounding boxes of the slices
rule create_bbox_extracts:
    input:
        "{OUTPUT_DIR}/json/{DATASET}.json",
    output:
        "{OUTPUT_DIR}/json/{DATASET}_extracts.geojson",
    script:
        "../scripts/prepare-extracts.py"


"""
Test with:
snakemake --cores all results/json/tanzania-mini_extracts.geojson
"""
