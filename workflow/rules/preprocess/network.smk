"""
Process gridfinder elements for each box
"""

rule create_power_network:
    conda: "../../../environment.yml"
    input:
        plants="{OUTPUT_DIR}/power/slice/{BOX}/network/powerplants.geoparquet",
        targets="{OUTPUT_DIR}/power/slice/{BOX}/network/targets.geoparquet",
        gridfinder="{OUTPUT_DIR}/power/slice/{BOX}/network/gridfinder.geoparquet",
    output:
        edges="{OUTPUT_DIR}/power/slice/{BOX}/network/edges.geoparquet",
        nodes="{OUTPUT_DIR}/power/slice/{BOX}/network/nodes.geoparquet",
    script:
        "../../scripts/process/process_power_4_network.py"
