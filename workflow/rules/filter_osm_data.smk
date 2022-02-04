# Take a .osm.pbf file and return a .osm.pbf file with a subset of the information
rule filter_osm_data:
    input:
        config["osmium_tags_filters_file"],
        os.path.join(config['data_dir'], f"{config['dataset']}.osm.pbf"),
    output:
        os.path.join(f"{config['output_dir']}",f"{config['dataset']}.highway-core.osm.pbf"),
    shell:
        "osmium tags-filter {input[1]} w/highway=$(cat {input[0]}) -o {output}"