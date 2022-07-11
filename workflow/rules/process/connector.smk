"""Finds connector points for all boxes

"""

out_connector = (
    expand(
        os.path.join(
            config["output_dir"],
            "power_processed",
            "all_boxes",
            "{box_id}",
            "connector_{box_id}.txt",
        ),
        box_id=all_boxes,
    ),
)


rule process_connector:
    input:
        expand(
            os.path.join(
                config["output_dir"],
                "power_processed",
                "all_boxes",
                "{box_id}",
                "network_{box_id}.gpkg",
            ),
            box_id=all_boxes,
        ),
        os.path.join(
            config["output_dir"], "power_processed", "world_boxes_metadata.txt"
        ),
    params:
        output_dir=config["output_dir"],
    output:
        out_connector,
    script:
        os.path.join("..", "..", "scripts", "process", "process_power_5_connector.py")
