# Create *-extracts.json file that determines the bounding boxes of the slices
rule create_slice_rules:
    input:
        os.path.join(f"{config['output_dir']}", "json", f"{config['dataset']}.json"),
    output:
        os.path.join(f"{config['output_dir']}", "json", f"{config['dataset']}-extracts.geojson"),
    script:
        "../scripts/prepare-extracts.py"
