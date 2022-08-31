"""Split the world into boxes

"""


all_box_geoms = expand(
    os.path.join(
        config["output_dir"],
        "power_processed",
        "all_boxes",
        "{box_id}",
        "geom_{box_id}.gpkg",
    ),
    box_id=ALL_BOXES,
)


rule world_splitter:
    input:
        os.path.join(
            config["output_dir"], "input", "admin-boundaries", "gadm36_levels.gpkg"
        ),
    output:
        all_box_geoms,
        os.path.join(
            config["output_dir"], "power_processed", "world_boxes_metadata.txt"
        ),
        os.path.join(config["output_dir"], "power_processed", "world_boxes.gpkg"),
    params:
        boxlen_value=config["box_width_height"],
        output_dir=config["output_dir"],
    script:
        os.path.join("..", "..", "scripts", "process", "world_split.py")
