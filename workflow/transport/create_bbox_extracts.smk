# Create *-extracts.json file that determines the bounding boxes of the slices
rule create_bbox_extracts:
    input:
        "{OUTPUT_DIR}/json/{DATASET}.json",
    params:
        # include slice_count as a param (despite not using elsewhere in the
        # rule) to trigger re-runs on change to this configuration option
        slice_count = config["slice_count"]
    output:
        # double curly braces allow us to expand but keep wildcards!
        expand(
            "{{OUTPUT_DIR}}/json/{{DATASET}}_extracts/slice-{n}.geojson",
            n=range(config["slice_count"]),
        ),
    script:
        "./prepare-extracts.py"


"""
Test with:
snakemake --cores all results/json/tanzania-mini_extracts.geojson
"""
