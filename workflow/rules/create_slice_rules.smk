# Create *-extracts.json file that determines the bounding boxes of the slices
rule create_slice_rules:
    input:
        os.path.join(f"{config['output_dir']}", "json", f"{dataset_slug}.json"),
    output:
        os.path.join(f"{config['output_dir']}", "json", f"{dataset_slug}_filter-{filter_slug}-extracts.geojson"),
    script:
        "../scripts/prepare-extracts.py"

rule test_create_slice_rules:
    input:
        os.path.join(f"{config['output_dir']}", "json", f"{dataset_slug}-{filter_slug}-extracts.geojson")