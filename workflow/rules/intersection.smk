# Combine .geoparquet files with .tiff files and output .geoparquet and .parquet files

from glob import glob

rule intersection:
    input:
        network="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}_{SLICE_SLUG}.geoparquet",
        # We read in the entire directory here to avoid splitting the job for each *.tif file
        tifs= lambda wildcards: checkpoints.trim_hazard_data.get(**wildcards).output[0],
    output:
        geoparquet="{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}_{SLICE_SLUG}_{HAZARD_SLUG}.geoparquet",
        parquet="{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}_{SLICE_SLUG}_{HAZARD_SLUG}.parquet",
    script:
        "../scripts/intersection.py"

"""
Test with:
snakemake --cores all results/splits/wales-latest_filter-highway-core_slice-0_hazard-aqueduct-coast.geoparquet
"""
