rule convert_to_geoparquet:
    input:
        os.path.join(DATA_DIR, "slices", "{slug}.highway-core.osm.pbf")
    output:
        os.path.join(DATA_DIR, "slices", "{slug}.highway-core.geoparquet")
    script:
        "osm_to_pq.py"
