# Take .geoparquet files of network and output .geoparquet files of network annotated with flow data

rule annotate_road_network:
    input:
        nodes="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}_road_nodes.geoparquet",
        edges="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}_road_edges.geoparquet",
        admin="{OUTPUT_DIR}/input/admin-boundaries/gadm36_levels.gpkg",
    output:
        nodes="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}_road_nodes_annotated.geoparquet",
        edges="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}_road_edges_annotated.geoparquet",
    script:
        "../scripts/annotate_road_network.py"

"""
Test with:
snakemake --cores all results/geoparquet/tanzania-mini_filter-highway-core/slice-0_road_edges_annotated.geoparquet
"""
