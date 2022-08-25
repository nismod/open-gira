# Use osmium to cut a .osm.pbf file into several .osm.pbf files as defined by bounding boxes


# This is a checkpoint because the DAG needs to be recalculated after the files are produced
checkpoint slice:
    input:
        data="{OUTPUT_DIR}/input/{DATASET}_{FILTER_SLUG}.osm.pbf",
        extracts_config="{OUTPUT_DIR}/json/{DATASET}_extracts.geojson",
    output:
        directory("{OUTPUT_DIR}/slices/{DATASET}_{FILTER_SLUG}"),
    run:
        import json
        import os
        from subprocess import run

        with open(input.extracts_config, "r") as fp:
            conf = json.load(fp)

        run(["mkdir", "-p", output[0]])
        for extract in conf["extracts"]:
            run([
                "osmium",
                "extract",
                "--set-bounds",
                "-b",
                ",".join([str(coord) for coord in extract["bbox"]]),
                input.data,
                "-o",
                os.path.join(output[0], extract["output"])
            ])




"""
Test with:
snakemake --cores all results/slices/tanzania-mini_filter-highway-core

Note: even testing with a single slice will generate all slices defined in the *_extracts.geojson file
"""
