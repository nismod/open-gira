# Create .json file that determines the bounding boxes of the complete data
rule create_overall_bbox:
    input:
        "{OUTPUT_DIR}/input/{DATASET}.osm.pbf",
        # os.path.join(f"{config['output_dir']}", "input", f"{dataset_slug}_filter-{filter_slug}.osm.pbf"),
    output:
        "{OUTPUT_DIR}/json/{DATASET}.json",
        # os.path.join(f"{config['output_dir']}", "json", f"{dataset_slug}.json"),
    script:
        "../../scripts/create_overall_bbox.py"


"""
Test with:
snakemake --cores all results/json/tanzania-mini.json
"""
