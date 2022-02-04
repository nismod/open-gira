# Use osmium to cut a .osm.pbf file into several .osm.pbf files as defined by bounding boxes
rule slice:
    input:
        # data=os.path.join(DATA_DIR, f"{DATASET}.osm.pbf"),
        data=os.path.join(f"{config['output_dir']}",f"{config['dataset']}.highway-core.osm.pbf"),
        extracts_config=os.path.join(f"{config['output_dir']}", "json", f"{config['dataset']}-extracts.geojson"),
    output:
        expand(
            os.path.join(f"{config['output_dir']}", "slices", f"{config['dataset']}-slice{{i}}.highway-core.osm.pbf"),
            i=range(config['slice_count'])
        ),
    shell:
        """
        osmium extract --overwrite --no-progress --config {input.extracts_config} {input.data}
        """
