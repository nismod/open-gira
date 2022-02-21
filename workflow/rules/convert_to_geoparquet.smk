# Take .osm.pbf files and output .geoparquet files

# https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#data-dependent-conditional-execution
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.slice.get(**wildcards).output[0]
    # Handle the slice wildcard manually, because it's not set by the slice.smk rule
    input = expand(os.path.join(checkpoint_output, "{SLICE_SLUG}.osm.pbf"), SLICE_SLUG=wildcards.SLICE_SLUG)
    print(f"checkpoint_output={checkpoint_output}")
    print(f"input={input}")
    return input

rule convert_to_geoparquet:
    input:
        aggregate_input,
    output:
        "{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}.geoparquet",
    script:
        "../scripts/osm_to_pq.py"

"""
Test with:
snakemake --cores all results/geoparquet/tanzania-mini_filter-highway-core/slice-0.geoparquet
"""