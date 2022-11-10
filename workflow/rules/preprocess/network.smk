"""
Process gridfinder elements for each box
"""

rule process_network:
    conda: "../../../environment.yml"
    input:
        plants="{OUTPUT_DIR}/processed/power/{BOX}/powerplants_{BOX}.parquet",
        targets="{OUTPUT_DIR}/processed/power/{BOX}/targets_{BOX}.parquet",
        gridfinder="{OUTPUT_DIR}/processed/power/{BOX}/gridfinder_{BOX}.parquet",
    output:
        edges="{OUTPUT_DIR}/processed/power/{BOX}/edges_{BOX}.parquet",
        nodes="{OUTPUT_DIR}/processed/power/{BOX}/nodes_{BOX}.parquet",
    script:
        "../../scripts/process/process_power_4_network.py"
