"""
Process gridfinder elements for each box
"""


rule process_gridfinder:
    input:
        os.path.join(config["output_dir"], "input", "gridfinder", "grid.gpkg"),
        os.path.join(
            config["output_dir"], "power_processed", "world_boxes_metadata.json"
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
                "gridfinder_{box_id}.gpkg",
            ),
            box_id=ALL_BOXES,
        )
    script:
        os.path.join("..", "..", "scripts", "process", "process_power_3_gridfinder.py")
