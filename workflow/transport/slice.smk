# Use osmium to cut a .osm.pbf file into several .osm.pbf files as defined by bounding boxes
rule slice:
    input:
        data="{OUTPUT_DIR}/input/OSM/{DATASET}_{FILTER_SLUG}.osm.pbf",
        extracts_config="{OUTPUT_DIR}/json/{DATASET}_extracts/{SLICE_SLUG}.geojson",
    output:
        slice_path="{OUTPUT_DIR}/slices/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}.osm.pbf",
    resources:
        mem_mb=32_000
    shell:
        # Need to run osmium from different output folders depending on the
        # FILTER_SLUG, and the extracts_config specifies only the filename.
        """
        WORKDIR=$(dirname {output.slice_path})
        mkdir -p $WORKDIR &&
        cd $WORKDIR &&
        osmium extract --set-bounds --overwrite --no-progress --config ../../../{input.extracts_config} ../../../{input.data}
        """


"""
Test with:
snakemake --cores all results/slices/tanzania-mini_filter-road/slice-0.osm.pbf
"""
