rule filter_osm_data:
    input:
        config["osmium_tags_filters_file"],
        os.path.join(DATA_DIR, "slices", "{slug}.osm.pbf")
    output:
        os.path.join(DATA_DIR, "slices", "{slug}.highway-core.osm.pbf")
    shell:
        "osmium tags-filter {input[1]} w/highway=$(cat {input[0]}) -o {output}"
