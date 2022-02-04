# Combine .geoparquet files with .tiff files and output .geoparquet and .parquet files
rule network_hazard_intersection:
    input:
        network="{OUTPUT_DIR}/geoparquet/{SLUG}.highway-core.geoparquet",
        tifs=glob(os.path.join(config['hazard_data_dir'], "*.tif")),
        hazard_csv=config["hazard_csv"],
    output:
        geoparquet="{OUTPUT_DIR}/splits/{SLUG}.highway-core_{HAZARD_SLUG}_splits.geoparquet",
        parquet="{OUTPUT_DIR}/splits/{SLUG}.highway-core_{HAZARD_SLUG}_splits.parquet",
    script:
        "../scripts/network_hazard_intersection.py"
