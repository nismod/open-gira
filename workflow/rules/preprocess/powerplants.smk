"""
Process powerplants for each box
"""


rule process_powerplants:
    conda: "../../../environment.yml"
    input:
        powerplants="{OUTPUT_DIR}/input/powerplants/global_power_plant_database.csv",
    output:
        powerplants="{OUTPUT_DIR}/processed/powerplants.geoparquet",
    script:
        "../../scripts/process/process_power_1_powerplants.py"
