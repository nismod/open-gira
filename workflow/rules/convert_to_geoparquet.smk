# Take .osm.pbf files and output .geoparquet files

# https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#data-dependent-conditional-execution
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.slice.get(**wildcards).output[0]
    # Handle the slice wildcard manually, because it's not set by the slice.smk rule
    return os.path.join(checkpoint_output, f"{wildcards.SLICE_SLUG}.osm.pbf")


rule convert_to_geoparquet:
    input:
        pbf=aggregate_input,
    params:
        # use the FILTER_SLUG to lookup the relevant keep_tags for this network type
        keep_tags=lambda wildcards: config["keep_tags"][wildcards.FILTER_SLUG.replace('filter-', '')]
    output:
        edges="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/raw/{SLICE_SLUG}_edges.geoparquet",
        nodes="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/raw/{SLICE_SLUG}_nodes.geoparquet",
    script:
        "../scripts/osm_to_pq.py"


"""
Test with:
snakemake --cores all results/geoparquet/tanzania-mini_filter-highway-core/raw/slice-0_edges.geoparquet
"""
