# Take .geoparquet files and output .geopackage files of annotated networks

rule clean_and_annotate_road_network:
    input:
        "{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}.geoparquet",
    output:
        "{OUTPUT_DIR}/geopackage/{DATASET}_{FILTER_SLUG}/{SLICE_SLUG}.gpkg",
    script:
        "../scripts/clean_and_annotate_pq.py"

"""
Test with:
snakemake --cores all results/geopackage/tanzania-mini_filter-highway-core/slice-0.gpkg
"""
