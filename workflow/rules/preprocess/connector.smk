"""
Finds connector points for all boxes
"""

rule process_connector:
    conda: "../../../environment.yml"
    input:
        edges="{OUTPUT_DIR}/processed/power/{BOX}/edges_{BOX}.parquet",
        nodes="{OUTPUT_DIR}/processed/power/{BOX}/nodes_{BOX}.parquet",
        global_metadata=rules.world_splitter.output.global_metadata,
    output:
        connector="{OUTPUT_DIR}/processed/power/{BOX}/connector_{BOX}.json",
    script:
        "../../scripts/process/process_power_5_connector.py"
