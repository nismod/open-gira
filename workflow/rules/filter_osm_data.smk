# Take a .osm.pbf file and return a .osm.pbf file with a subset of the information
rule filter_osm_data:
    input:
        config["osmium_tags_filters_file"],
        os.path.join(config['data_dir'], f"{config['dataset']}.osm.pbf"),
    output:
        "{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}.osm.pbf",
        # os.path.join(f"{config['output_dir']}",f"{config['dataset']}_filter-{filter_slug}.osm.pbf"),
    shell:
        "osmium tags-filter {input[1]} w/highway=$(cat {input[0]}) -o {output}"

rule test_filter_osm_data:
    input:
        os.path.join(f"{config['output_dir']}",f"{config['dataset']}_filter-{filter_slug}.osm.pbf")