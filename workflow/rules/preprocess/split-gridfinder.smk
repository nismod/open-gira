"""
Process gridfinder elements for each box
"""


rule process_gridfinder:
    input:
        os.path.join(config["output_dir"], "input", "gridfinder", "grid.gpkg"),
        rules.world_splitter.output.global_metadata,
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
