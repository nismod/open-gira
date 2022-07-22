# Take .geoparquet files of network and output .geoparquet files of network annotated with flow data


rule annotate_transport_network:
    input:
        nodes="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/network/{SLICE_SLUG}_nodes.geoparquet",
        edges="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/network/{SLICE_SLUG}_edges.geoparquet",
        admin="{OUTPUT_DIR}/input/admin-boundaries/gadm36_levels.gpkg",
    output:
        nodes="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/network/{SLICE_SLUG}_nodes_annotated.geoparquet",
        edges="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/network/{SLICE_SLUG}_edges_annotated.geoparquet",
    params:
        # determine the network type from the filter, e.g. road, rail
        network_type=lambda wildcards: wildcards.FILTER_SLUG.replace('filter-', '')
    script:
        # template the path string with a value from params (can't execute .replace in `script` context)
        "../scripts/annotate_{params.network_type}_network.py"


"""
Test with:
snakemake --cores all results/geoparquet/tanzania-mini_filter-road/slice-0_edges_annotated.geoparquet
"""
