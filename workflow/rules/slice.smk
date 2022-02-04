rule slice:
    input:
        # data=os.path.join(DATA_DIR, f"{DATASET}.osm.pbf"),
        data=os.path.join(f"{OUTPUT_DIR}",f"{DATASET}.highway-core.osm.pbf"),
        extracts_config=os.path.join(f"{OUTPUT_DIR}", "json", f"{DATASET}-extracts.geojson"),
    output:
        expand(
            os.path.join(f"{OUTPUT_DIR}", "slices", f"{DATASET}-slice{{i}}.highway-core.osm.pbf"),
            i=range(config['slice_count'])
        ),
    shell:
        """
        osmium extract --overwrite --no-progress --config {input.extracts_config} {input.data}
        """
