rule parse_IRIS:
    input:
        csv_dir="{OUTPUT_DIR}/input/IRIS/iris-data/event_sets/{IRIS_SCENARIO}/"
    output:
        parquet="{OUTPUT_DIR}/storm_tracks/IRIS-{IRIS_SCENARIO}/{SAMPLE}/tracks.geoparquet"
    script:
        "../../scripts/preprocess/parse_IRIS.py"

"""
Test with:
snakemake -c1 results/storm_tracks/IRIS_SSP1-2050/0/tracks.geoparquet
"""


rule slice_IRIS:
    input:
        global_tracks=rules.parse_IRIS.output.parquet,
        grid_hull="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/convex_hull.json"
    output:
        sliced_tracks="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/IRIS-{IRIS_SCENARIO}/{SAMPLE}/tracks.geoparquet",
    resources:
        mem_mb=10000  # the global tracks file is fairly chunky
    script:
        "../../scripts/preprocess/slice_storm_tracks.py"

"""
To test:
snakemake -c1 results/power/by_country/PRI/storms/IRIS-PRESENT/0/tracks.geoparquet
"""
