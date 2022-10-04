"""
Finds connector points for all boxes
"""

rule process_connector:
    conda: "../../../environment.yml"
    input:
        expand(
            os.path.join(
                config["output_dir"],
                "power_processed",
                "all_boxes",
                "{box_id}",
                "network_{box_id}.gpkg",
            ),
            box_id=ALL_BOXES,
        ),
        os.path.join(
            config["output_dir"], "power_processed", "world_boxes_metadata.txt"
        ),
    params:
        output_dir=config["output_dir"],
    output:
        CONNECTOR_OUT,
    script:
        os.path.join("..", "..", "scripts", "process", "process_power_5_connector.py")
