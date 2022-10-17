# Combine .geoparquet files with .tiff files and output .geoparquet and .parquet files

from glob import glob


rule network_raster:
    conda: "../../../environment.yml"
    input:
        network="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/processed/{SLICE_SLUG}_edges.geoparquet",
        # We read in the entire directory here to avoid splitting the job for each *.tif file
        tifs=lambda wildcards: checkpoints.trim_hazard_data.get(**wildcards).output.trimmed_dir,
    output:
        geoparquet="{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/{SLICE_SLUG}.geoparquet",
        parquet="{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/{SLICE_SLUG}.parquet",
    # intersection.py is intermittently failing, use retries as temporary fix
    # TODO: investigate why this job is sometimes failing with a coredump from intersection.py
    retries: 3
    script:
        "../../scripts/intersection.py"


"""
Test with:
snakemake --cores all results/splits/tanzania-mini_filter-road/hazard-aqueduct-river/slice-0.geoparquet
"""
