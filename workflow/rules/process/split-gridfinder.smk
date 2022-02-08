"""Process gridfinder elements for each box

"""

out_gridfinder = (
    expand(
        os.path.join(
            DATA_DIR, "processed", "all_boxes", "{box_id}", "gridfinder_{box_id}.gpkg"
        ),
        box_id=all_boxes,
    ),
)


rule process_gridfinder:
    input:
        os.path.join("data", "gridfinder", "grid.gpkg"),
        os.path.join("data", "processed", "world_boxes_metadata.txt"),
    output:
        out_gridfinder,
    shell:
        "python3 " + os.path.join(
        "workflow", "scripts", "process", "process_power_3_gridfinder.py"
        )
