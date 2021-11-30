rule convert_to_geoparquet:
    input:
        "results/filtered/{slug}.highway-core.osm.pbf"
    output:
        "results/geoparquet/{slug}.highway-core.geoparquet"
    script:
        "../scripts/osm_to_pq.py"
