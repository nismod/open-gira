"""Finds connector points for all boxes

"""

out_connector = (
    expand(
        os.path.join(
            "data", "processed", "all_boxes", "{box_id}", "connector_{box_id}.txt"
        ),
        box_id=all_boxes,
    ),
)


rule process_connector:
    input:
        expand(
            os.path.join(
                "data", "processed", "all_boxes", "{box_id}", "network_{box_id}.gpkg"
            ),
            box_id=all_boxes,
        ),
        os.path.join("data", "processed", "world_boxes_metadata.txt"),
    output:
        out_connector,
    shell:
        "python3 " + os.path.join(
        "workflow", "scripts", "process", "process_power_5_connector.py"
        )
