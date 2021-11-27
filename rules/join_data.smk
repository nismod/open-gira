def aggregate_input_geoparquet(wildcards):
    checkpoint_output = checkpoints.slice.get(**wildcards).output[0]
    return expand(
        os.path.join(
            OUTPUT_DIR,
            "slices",
            f"{DATASET}-slice{{i}}.highway-core_{hazard_slug}_splits.geoparquet",
        ),
        i=glob_wildcards(os.path.join(checkpoint_output, f"{DATASET}-slice{{i,\d+}}.osm.pbf")).i
    )

rule join_data:
    input:
        aggregate_input_geoparquet
    output:
        os.path.join(OUTPUT_DIR, f"{DATASET}.highway-core_{hazard_slug}_splits.geoparquet")
    script:
        "join_data.py"
