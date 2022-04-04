"""Process powerplants for each box

"""

out_powerplants = (
    expand(
        os.path.join(
            'data', "processed", "all_boxes", "{box_id}", "powerplants_{box_id}.csv"
        ),
        box_id=all_boxes,
    ),
)


rule process_powerplants:
    input:
        os.path.join('data', "powerplants", "global_power_plant_database.csv"),
        os.path.join("data", "processed", "world_boxes_metadata.txt"),
    output:
        out_powerplants,
    script:
            os.path.join("..", "..", "scripts", "process", "process_power_1_powerplants.py"
        )
