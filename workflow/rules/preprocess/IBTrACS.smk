rule parse_ibtracs:
    input:
        ibtracs_csv = "results/input/IBTrACS/raw/v4.csv"
    output:
        ibtracs_parquet = "results/input/IBTrACS/processed/v4.geoparquet"
    script:
        "../../scripts/process/parse_IBTrACS.py"

"""
To test:
snakemake -c1 results/input/IBTrACS/processed/v4.geoparquet
"""


rule slice_ibtracs:
    input:
        global_tracks=rules.parse_ibtracs.output.ibtracs_parquet,
        global_boxes=rules.world_splitter.output.global_boxes,
    output:
        sliced_tracks="{OUTPUT_DIR}/power/slice/{BOX}/storms/IBTrACS.geoparquet",
    script:
        "../../scripts/process/slice_storm_tracks.py"

"""
To test:
snakemake -c1 results/power/slice/1030/storms/IBTrACS.geoparquet
"""
