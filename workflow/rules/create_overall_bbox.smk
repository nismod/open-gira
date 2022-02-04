# Create .json file that determines the bounding boxes of the complete data
rule create_overall_bbox:
    input:
        os.path.join(f"{config['output_dir']}", f"{config['dataset']}.highway-core.osm.pbf"),
    output:
        os.path.join(f"{config['output_dir']}", "json", f"{config['dataset']}.json"),
    script:
        "../scripts/create_overall_bbox.py"
