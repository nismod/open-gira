# Take .geoparquet OSM files and output files of cleaned network nodes and edges

rule create_transport_network:
    input:
        "{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}.geoparquet"
    output:
        # can we parameterise the network type (i.e. road) in the output filename?
        # config["transport_type"] is presumably available
        "{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}_road_nodes.geoparquet",
        "{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}_road_edges.geoparquet"
    script:
        "../scripts/create_transport_network.py"

"""
Test with:
snakemake --cores all results/geoparquet/tanzania-mini_filter-highway-core/slice-0_road_edges.geoparquet
"""
