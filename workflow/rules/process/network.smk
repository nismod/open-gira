"""Process gridfinder elements for each box

"""

out_network = (
    expand(
        [
            os.path.join(
                config["output_dir"],
                "power_processed",
                "all_boxes",
                "{box_id}",
                "plants_{box_id}.gpkg",
            ),
            os.path.join(
                config["output_dir"],
                "power_processed",
                "all_boxes",
                "{box_id}",
                "network_{box_id}.gpkg",
            ),
            os.path.join(
                config["output_dir"],
                "power_processed",
                "all_boxes",
                "{box_id}",
                "targets_{box_id}.gpkg",
            ),
        ],
        box_id=all_boxes,
    ),
)


rule process_network:
    input:
        os.path.join(
            config["output_dir"],
            "power_processed",
            "all_boxes",
            "{box_id}",
            "powerplants_{box_id}.csv",
        ),
        os.path.join(
            config["output_dir"],
            "power_processed",
            "all_boxes",
            "{box_id}",
            "targets_{box_id}.csv",
        ),
        os.path.join(
            config["output_dir"],
            "power_processed",
            "all_boxes",
            "{box_id}",
            "gridfinder_{box_id}.gpkg",
        ),
        os.path.join(
            config["output_dir"], "power_processed", "world_boxes_metadata.txt"
        ),
    output:
        os.path.join(
            config["output_dir"],
            "power_processed",
            "all_boxes",
            "{box_id}",
            "plants_{box_id}.gpkg",
        ),
        os.path.join(
            config["output_dir"],
            "power_processed",
            "all_boxes",
            "{box_id}",
            "network_{box_id}.gpkg",
        ),
        os.path.join(
            config["output_dir"],
            "power_processed",
            "all_boxes",
            "{box_id}",
            "targets_{box_id}.gpkg",
        ),
    params:
        box_id="{box_id}",
        output_dir=config["output_dir"],
    script:
        os.path.join("..", "..", "scripts", "process", "process_power_4_network.py")
