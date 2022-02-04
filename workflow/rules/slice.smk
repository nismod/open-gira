rule slice:
    input:
        # data=os.path.join(DATA_DIR, f"{DATASET}.osm.pbf"),
        data=os.path.join(f"{OUTPUT_DIR}",f"{DATASET}.highway-core.osm.pbf"),
        extracts_config=os.path.join(f"{DATA_DIR}", f"{DATASET}-extracts.geojson"),
    output:
        slices=directory(f"{OUTPUT_DIR}/slices"),
    shell:
        """
        mkdir {output.slices} &&
        osmium extract --no-progress --config {input.extracts_config} {input.data} \
            --directory {output.slices}
        """
