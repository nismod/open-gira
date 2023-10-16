rule parse_ibtracs:
    input:
        ibtracs_csv = "results/input/IBTrACS/raw/v4.csv"
    output:
        ibtracs_parquet = "results/storm_tracks/IBTrACS/0/tracks.geoparquet"
    script:
        "../../scripts/preprocess/parse_IBTrACS.py"

"""
To test:
snakemake -c1 results/storm_tracks/IBTrACS/0/tracks.geoparquet
"""


rule slice_ibtracs:
    input:
        global_tracks=rules.parse_ibtracs.output.ibtracs_parquet,
        grid_hull="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/convex_hull.json"
    output:
        sliced_tracks="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/IBTrACS/0/tracks.geoparquet",
    script:
        "../../scripts/preprocess/slice_storm_tracks.py"

"""
To test:
snakemake -c1 results/power/by_country/PRI/storms/IBTrACS/0/tracks.geoparquet
"""
