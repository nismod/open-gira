# Take .osm.pbf files and output .geoparquet files
rule convert_to_geoparquet:
    input:
        "{OUTPUT_DIR}/slices/{SLUG}.osm.pbf",
    output:
        "{OUTPUT_DIR}/geoparquet/{SLUG}.geoparquet",
    script:
        "../scripts/osm_to_pq.py"
