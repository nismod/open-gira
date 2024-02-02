# Take .osm.pbf files and output .geoparquet files
rule osm_to_geoparquet:
    input:
        pbf=rules.slice.output.slice_path,
    params:
        # use the FILTER_SLUG to lookup the relevant keep_tags for this network type
        keep_tags=lambda wildcards: config["keep_tags"][wildcards.FILTER_SLUG.replace('filter-', '')]
    output:
        edges="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/raw/{SLICE_SLUG}_edges.geoparquet",
        nodes="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/raw/{SLICE_SLUG}_nodes.geoparquet",
    script:
        "./osm_to_pq.py"


"""
Test with:
snakemake --cores all results/geoparquet/tanzania-mini_filter-road/raw/slice-0_edges.geoparquet
"""
