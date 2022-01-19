def aggregate_input_geoparquet(wildcards):
    checkpoint_output = checkpoints.slice.get(**wildcards).output[0]
    return expand(
        f"{OUTPUT_DIR}/splits/{DATASET}-slice{{i}}.highway-core_{hazard_slug}_splits.geoparquet",
        i=glob_wildcards(
            os.path.join(checkpoint_output, f"{DATASET}-slice{{i,\d+}}.osm.pbf")
        ).i,
    )


rule join_data:
    input:
        aggregate_input_geoparquet,
    output:
        f"{OUTPUT_DIR}/{DATASET}.highway-core_{hazard_slug}_splits.geoparquet",
    script:
        "../scripts/join_data.py"
