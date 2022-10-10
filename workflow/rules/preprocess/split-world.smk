"""
Split the world into boxes
"""


rule world_splitter:
    conda: "../../../environment.yml"
    input:
        admin_data = rules.download_gadm_levels.output.admin_bounds_global_layer_per_level,
    output:
        geometry_by_box = expand(
            os.path.join(
                config["output_dir"],
                "power_processed",
                "all_boxes",
                "{box_id}",
                "geom_{box_id}.gpkg",
            ),
            box_id=ALL_BOXES,
        ),
        global_metadata = os.path.join(
            config["output_dir"], "power_processed", "world_boxes_metadata.json"
        ),
        global_boxes = os.path.join(config["output_dir"], "power_processed", "world_boxes.gpkg"),
    params:
        boxlen_value = config["box_width_height"],
        output_dir = config["output_dir"],
    script:
        os.path.join("..", "..", "scripts", "process", "world_split.py")
