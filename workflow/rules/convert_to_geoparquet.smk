rule convert_to_geoparquet:
    input:
        "{OUTPUT_DIR}/filtered/{slug}.highway-core.osm.pbf",
    output:
        "{OUTPUT_DIR}/geoparquet/{slug}.highway-core.geoparquet",
    script:
        "../scripts/osm_to_pq.py"
