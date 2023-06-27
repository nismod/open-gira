rule parse_storm:
    input:
        csv_dir="{OUTPUT_DIR}/input/STORM/events/{STORM_MODEL}/raw"
    output:
        parquet="{OUTPUT_DIR}/storm_tracks/STORM-{STORM_MODEL}/tracks.geoparquet"
    script:
        "../../scripts/preprocess/parse_STORM.py"

"""
Test with:
snakemake -c1 results/storm_tracks/STORM-constant/tracks.geoparquet
"""


rule slice_storm:
    input:
        global_tracks=rules.parse_storm.output.parquet,
        grid_hull="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/convex_hull.json"
    output:
        sliced_tracks="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/STORM-{STORM_MODEL}/tracks.geoparquet",
    resources:
        mem_mb=10000  # the global tracks file is fairly chunky
    script:
        "../../scripts/preprocess/slice_storm_tracks.py"

"""
To test:
snakemake -c1 results/power/by_country/PRI/storms/STORM-constant/tracks.geoparquet
"""
