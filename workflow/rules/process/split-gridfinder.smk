"""Process gridfinder elements for each box

"""

out_gridfinder = (
    expand(
        os.path.join(
            config['output_dir'], "power_processed", "all_boxes", "{box_id}", "gridfinder_{box_id}.gpkg"
        ),
        box_id=all_boxes,
    ),
)


rule process_gridfinder:
    input:
        os.path.join(config['output_dir'], "input", "gridfinder", "grid.gpkg"),
        os.path.join(config['output_dir'], "power_processed", "world_boxes_metadata.txt"),
    params:
        output_dir = config['output_dir']
    output:
        out_gridfinder,
    script:
        os.path.join("..", "..", "scripts", "process", "process_power_3_gridfinder.py")
