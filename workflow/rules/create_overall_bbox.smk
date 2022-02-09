# Create .json file that determines the bounding boxes of the complete data
rule create_overall_bbox:
    input:
        os.path.join(f"{config['output_dir']}", "input", f"{dataset_slug}_filter-{filter_slug}.osm.pbf"),
    output:
        os.path.join(f"{config['output_dir']}", "json", f"{dataset_slug}.json"),
    script:
        "../scripts/create_overall_bbox.py"

rule test_create_overall_bbox:
    input:
        os.path.join(f"{config['output_dir']}","json",f"{dataset_slug}.json")