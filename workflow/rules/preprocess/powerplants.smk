"""
Process powerplants for each box
"""


rule process_powerplants:
    conda: "../../../environment.yml"
    input:
        os.path.join(
            config["output_dir"],
            "input",
            "powerplants",
            "global_power_plant_database.csv",
        ),
        os.path.join(
            config["output_dir"], "power_processed", "world_boxes_metadata.txt"
        ),
    params:
        output_dir=config["output_dir"],
    output:
        expand(
            os.path.join(
                config["output_dir"],
                "power_processed",
                "all_boxes",
                "{box_id}",
                "powerplants_{box_id}.csv",
            ),
            box_id=ALL_BOXES,
        )
    script:
        os.path.join("..", "..", "scripts", "process", "process_power_1_powerplants.py")
