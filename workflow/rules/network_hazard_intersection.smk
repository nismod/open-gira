# Combine .geoparquet files with .tiff files and output .geoparquet and .parquet files
rule network_hazard_intersection:
    input:
        network="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}_{SLICE_SLUG}.geoparquet",
        tifs=glob(os.path.join(f"{config['output_dir']}", "input", "hazard-*", "*.tif"))
    output:
        geoparquet="{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}_{SLICE_SLUG}_{HAZARD_SLUG}.geoparquet",
        parquet="{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}_{SLICE_SLUG}_{HAZARD_SLUG}.parquet",
    script:
        "../scripts/network_hazard_intersection.py"

# Test via target file:
# snakemake --cores all results/splits/tanzania-latest_filter-highway-core_slice-0_hazard-aqueduct-coast.geoparquet
