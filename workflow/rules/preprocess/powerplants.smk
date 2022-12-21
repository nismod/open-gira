"""
Process powerplants for each box
"""


rule process_powerplants:
    conda: "../../../environment.yml"
    input:
        powerplants="{OUTPUT_DIR}/input/powerplants/global_power_plant_database.csv",
        global_boxes=rules.world_splitter.output.global_boxes,
    output:
        powerplants="{OUTPUT_DIR}/power/slice/{BOX}/network/powerplants.geoparquet",
    script:
        "../../scripts/process/process_power_1_powerplants.py"
