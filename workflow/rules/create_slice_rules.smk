# Create *-extracts.json file that determines the bounding boxes of the slices
rule create_slice_rules:
    input:
        os.path.join(f"{OUTPUT_DIR}", "json", f"{DATASET}.json"),
    output:
        os.path.join(f"{OUTPUT_DIR}", "json", f"{DATASET}-extracts.geojson"),
    script:
        "../scripts/prepare-extracts.py"
