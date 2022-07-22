# Use osmium to cut a .osm.pbf file into several .osm.pbf files as defined by bounding boxes


# This is a checkpoint because the DAG needs to be recalculated after the files are produced
checkpoint slice:
    input:
        data="{OUTPUT_DIR}/input/{DATASET}_{FILTER_SLUG}.osm.pbf",
        extracts_config="{OUTPUT_DIR}/json/{DATASET}_extracts.geojson",
    output:
        directory("{OUTPUT_DIR}/slices/{DATASET}_{FILTER_SLUG}"),
    shell:
        # Need to run osmium from different output folders depending on the
        # FILTER_SLUG, and the extracts_config specifies only the filename.
        """
        mkdir {output} &&
        cd {output} &&
        osmium extract --set-bounds --overwrite --no-progress --config ../../../{input.extracts_config} ../../../{input.data}
        """


"""
Test with:
snakemake --cores all results/slices/tanzania-mini_filter-highway-core

Note: even testing with a single slice will generate all slices defined in the *_extracts.geojson file
"""
