# Create *-extracts.json file that determines the bounding boxes of the slices
rule create_bbox_extracts:
    conda: "../../../environment.yml"
    input:
        "{OUTPUT_DIR}/json/{DATASET}.json",
    output:
        # double curly braces allow us to expand but keep wildcards!
        expand(
            "{{OUTPUT_DIR}}/json/{{DATASET}}_extracts/slice-{n}.json",
            n=range(config["slice_count"]),
        ),
    script:
        "../../scripts/prepare-extracts.py"


"""
Test with:
snakemake --cores all results/json/tanzania-mini_extracts.geojson
"""
