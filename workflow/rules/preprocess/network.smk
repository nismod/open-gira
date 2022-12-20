"""
Process gridfinder elements for each box
"""

rule create_power_network:
    conda: "../../../environment.yml"
    input:
        plants="{OUTPUT_DIR}/power/slice/{BOX}/network/powerplants_{BOX}.geoparquet",
        targets="{OUTPUT_DIR}/power/slice/{BOX}/network/targets_{BOX}.geoparquet",
        gridfinder="{OUTPUT_DIR}/power/slice/{BOX}/network/gridfinder_{BOX}.geoparquet",
    output:
        edges="{OUTPUT_DIR}/power/slice/{BOX}/network/edges_{BOX}.geoparquet",
        nodes="{OUTPUT_DIR}/power/slice/{BOX}/network/nodes_{BOX}.geoparquet",
    script:
        "../../scripts/process/process_power_4_network.py"
