"""Process powerplants for each box

"""

out_powerplants = (
    expand(
        os.path.join(
            DATA_DIR, "processed", "all_boxes", "{box_id}", "powerplants_{box_id}.csv"
        ),
        box_id=all_boxes,
    ),
)


rule process_powerplants:
    input:
        os.path.join(DATA_DIR, "powerplants", "global_power_plant_database.csv"),
        os.path.join("data", "processed", "world_boxes_metadata.txt"),
    output:
        out_powerplants,
    shell:
        "python3 " + os.path.join(
        "workflow", "scripts", "processing", "process_power_1_powerplants.py"
        )
