# Take .geoparquet OSM files and output .geopackage files of cleaned network

rule create_transport_network:
    input:
        "{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}.geoparquet"
    output:
        # can we parameterise the network type (i.e. road) in the output filename?
        # config["transport_type"] is presumably available
        "{OUTPUT_DIR}/geopackage/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}_road_network.gpkg"
    script:
        "../scripts/create_transport_network.py"

"""
Test with:
snakemake --cores all results/geopackage/tanzania-mini_filter-highway-core/slice-0_road_network.gpkg
"""
