rule network_hazard_intersection:
    input:
        network=os.path.join(DATA_DIR, "slices", "{slug}.highway-core.geoparquet"),
        tifs=glob(os.path.join(HAZARD_DATA_DIR, "*.tif")),
        hazard_csv=config["hazard_csv"]
    output:
        geoparquet=os.path.join(OUTPUT_DIR, "slices", "{slug}.highway-core_{hazard_slug}_splits.geoparquet"),
        parquet=os.path.join(OUTPUT_DIR, "slices", "{slug}.highway-core_{hazard_slug}_splits.parquet")
    script:
        "network_hazard_intersection.py"
