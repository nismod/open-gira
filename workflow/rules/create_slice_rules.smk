# Create *-extracts.json file that determines the bounding boxes of the slices
rule create_slice_rules:
    input:
        os.path.join(f"{config['output_dir']}", "json", f"{config['dataset']}.json"),
    output:
        os.path.join(f"{config['output_dir']}", "json", f"{config['dataset']}_filter-{filter_slug}-extracts.geojson"),
    script:
        "../scripts/prepare-extracts.py"

rule test_create_slice_rules:
    input:
        os.path.join(f"{config['output_dir']}", "json", f"{config['dataset']}-{filter_slug}-extracts.geojson")