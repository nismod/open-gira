# Create .json file that determines the bounding boxes of the complete data
rule create_overall_bbox:
    input:
        os.path.join(f"{OUTPUT_DIR}", f"{DATASET}.highway-core.osm.pbf"),
    output:
        os.path.join(f"{OUTPUT_DIR}", "json", f"{DATASET}.json"),
    script:
        "../scripts/create_overall_bbox.py"
