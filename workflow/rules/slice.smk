rule slice:
    input:
        # data=os.path.join(DATA_DIR, f"{DATASET}.osm.pbf"),
        data=os.path.join(f"{OUTPUT_DIR}",f"{DATASET}.highway-core.osm.pbf"),
        extracts_config=os.path.join(f"{DATA_DIR}", f"{DATASET}-extracts.geojson"),
    output:
        "{OUTPUT_DIR}/slices/{slug}.osm.pbf",
    shell:
        """
        osmium extract --overwrite --no-progress --config {input.extracts_config} {input.data}
        """
