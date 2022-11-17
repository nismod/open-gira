"""
Process gridfinder elements for each box
"""

rule process_network:
    conda: "../../../environment.yml"
    input:
        plants="{OUTPUT_DIR}/processed/power/{BOX}/powerplants_{BOX}.geoparquet",
        targets="{OUTPUT_DIR}/processed/power/{BOX}/targets_{BOX}.geoparquet",
        gridfinder="{OUTPUT_DIR}/processed/power/{BOX}/gridfinder_{BOX}.geoparquet",
    output:
        edges="{OUTPUT_DIR}/processed/power/{BOX}/edges_{BOX}.geoparquet",
        nodes="{OUTPUT_DIR}/processed/power/{BOX}/nodes_{BOX}.geoparquet",
    script:
        "../../scripts/process/process_power_4_network.py"
