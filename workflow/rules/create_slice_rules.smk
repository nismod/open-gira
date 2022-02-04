# Create *-extracts.json file that determines the bounding boxes of the slices
rule create_slice_rules:
    input:
        os.path.join(f"{DATA_DIR}", f"{DATASET}.json"),
    output:
        os.path.join(f"{DATA_DIR}", f"{DATASET}-extracts.geojson"),
    script:
        "../scripts/prepare-extracts.py"
