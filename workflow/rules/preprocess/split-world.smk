"""Split the world into boxes

"""


## Process Config Inputs ##
# Preliminary Snakemake Calculations
all_boxes = [
    f"box_{int(idx)}"
    for idx in range(
        0, int((180 - -180) * (90 - -90) / config["box_width_height"] ** 2)
    )
]  # all boxes with width and height boxlen


if config["specific_boxes"] != "None":
    print("Using specified boxes")
    all_boxes = [f"box_{num}" for num in config["specific_boxes"]]
if len(all_boxes) == 0:
    print("Specific boxes incorrectly specified")


all_box_geoms = expand(
    os.path.join(
        config["output_dir"],
        "power_processed",
        "all_boxes",
        "{box_id}",
        "geom_{box_id}.gpkg",
    ),
    box_id=all_boxes,
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
