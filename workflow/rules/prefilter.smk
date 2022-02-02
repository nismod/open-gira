rule filter:
    input:
        config["osmium_tags_filters_file"],
        os.path.join(DATA_DIR, f"{DATASET}.osm.pbf"),
    output:
        os.path.join(f"{OUTPUT_DIR}",f"{DATASET}.highway-core.osm.pbf"),
    shell:
        "osmium tags-filter {input[1]} w/highway=$(cat {input[0]}) -o {output}"