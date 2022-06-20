# Take .geoparquet files of network and output .geoparquet files of network annotated with flow data

rule annotate_road_network:
    input:
        "{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}_road_nodes.geoparquet",
        "{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}_road_edges.geoparquet"
    output:
        "{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}_road_nodes_annotated.geoparquet",
        "{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}_road_edges_annotated.geoparquet"
    script:
        "../scripts/annotate_road_network.py"

"""
Test with:
snakemake --cores all results/geoparquet/tanzania-mini_filter-highway-core/slice-0_road_edges_annotated.geoparquet
"""
