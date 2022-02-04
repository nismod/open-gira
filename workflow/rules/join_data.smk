rule join_data:
    input:
        expand(
            os.path.join(
                f"{OUTPUT_DIR}",
                "splits",
                f"{DATASET}-slice{{i}}.highway-core_{hazard_slug}_splits.geoparquet"
            ),
            i=range(config['slice_count'])
        ),
    output:
        f"{OUTPUT_DIR}/{DATASET}.highway-core_{hazard_slug}_splits.geoparquet",
    script:
        "../scripts/join_data.py"
