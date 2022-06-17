# Take .geopackage files of network and output .geopackage files of network annotated with flow data

rule annotate_road_network:
    input:
        "{OUTPUT_DIR}/geopackage/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}_road_network.gpkg"
    output:
        "{OUTPUT_DIR}/geopackage/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}_road_network_annotated.gpkg"
    script:
        "../scripts/annotate_road_network.py"

"""
Test with:
snakemake --cores all results/geopackage/tanzania-mini_filter-highway-core/slice-0_road_network_annotated.gpkg
"""
