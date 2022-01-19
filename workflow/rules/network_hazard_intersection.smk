rule network_hazard_intersection:
    input:
        network="{OUTPUT_DIR}/geoparquet/{slug}.highway-core.geoparquet",
        tifs=glob(os.path.join(HAZARD_DATA_DIR, "*.tif")),
        hazard_csv=config["hazard_csv"],
    output:
        geoparquet="{OUTPUT_DIR}/splits/{slug}.highway-core_{hazard_slug}_splits.geoparquet",
        parquet="{OUTPUT_DIR}/splits/{slug}.highway-core_{hazard_slug}_splits.parquet",
    script:
        "../scripts/network_hazard_intersection.py"
