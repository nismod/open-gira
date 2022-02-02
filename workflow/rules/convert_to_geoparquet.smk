rule convert_to_geoparquet:
    input:
        "{OUTPUT_DIR}/slices/{slug}.osm.pbf",
    output:
        "{OUTPUT_DIR}/geoparquet/{slug}.highway-core.geoparquet",
    script:
        "../scripts/osm_to_pq.py"
