rule parse_storm:
    input:
        csv_dir=directory("{OUTPUT_DIR}/input/STORM/events/{STORM_MODEL}/raw")
    output:
        parquet="{OUTPUT_DIR}/input/STORM/events/{STORM_MODEL}/processed.parquet"
    script:
        "../../scripts/process/parse_STORM.py"

"""
Test with:
snakemake -c1 results/input/STORM/events/constant/processed.parquet
"""


rule slice_storm:
    input:
        global_tracks=rules.parse_storm.output.parquet,
        global_boxes=rules.world_splitter.output.global_boxes,
    output:
        sliced_tracks="{OUTPUT_DIR}/power/slice/{BOX}/storms/{STORM_MODEL}/STORM.geoparquet",
    script:
        "../../scripts/process/slice_storm_tracks.py"

"""
To test:
snakemake -c1 results/power/slice/1030/storms/constant/STORM.geoparquet
"""
