# Combine .geoparquet files with .tiff files and output .geoparquet and .parquet files

from glob import glob


rule intersection:
    input:
        network="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/processed/{SLICE_SLUG}_edges.geoparquet",
        # We read in the entire directory here to avoid splitting the job for each *.tif file
        tifs=lambda wildcards: checkpoints.trim_hazard_data.get(**wildcards).output[0],
    output:
        geoparquet="{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/{SLICE_SLUG}.geoparquet",
        parquet="{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/{SLICE_SLUG}.parquet",
    script:
        "../scripts/intersection.py"


"""
Test with:
snakemake --cores all results/splits/tanzania-mini_filter-highway-core/hazard-aqueduct-river/slice-0_edges.geoparquet
"""
