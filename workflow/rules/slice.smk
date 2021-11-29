checkpoint slice:
    input:
        data=os.path.join(DATA_DIR, f"{DATASET}.osm.pbf"),
        extracts_config=os.path.join(DATA_DIR, f"{DATASET}-extracts.geojson")
    output:
        slices=directory(os.path.join(DATA_DIR, "slices"))
    shell:
        """mkdir {output.slices} &&
        osmium extract --no-progress --config {input.extracts_config} {input.data} \
        --directory {output.slices}
        """
