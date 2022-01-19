rule filter_osm_data:
    input:
        config["osmium_tags_filters_file"],
        "results/slices/{slug}.osm.pbf",
    output:
        "results/filtered/{slug}.highway-core.osm.pbf",
    shell:
        "osmium tags-filter {input[1]} w/highway=$(cat {input[0]}) -o {output}"
