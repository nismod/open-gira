"""Simplify each box

"""

rule process_simplify:
    input:
        os.path.join(
            "data", "processed", "all_boxes", "{box_id}", "network_{box_id}.gpkg"
        ),
        os.path.join(
            "data", "processed", "all_boxes", "{box_id}", "connector_{box_id}.txt"
        ),
        os.path.join(
            "data", "processed", "all_boxes", "{box_id}", "plants_{box_id}.gpkg"
        ),
        os.path.join(
            "data", "processed", "all_boxes", "{box_id}", "targets_{box_id}.gpkg"
        ),
    output:
        os.path.join(
            "data",
            "processed",
            "all_boxes",
            "{box_id}",
            "simple_network_{box_id}.gpkg"
        ),
        os.path.join(
            "data",
            "processed",
            "all_boxes",
            "{box_id}",
            "collapsed_sources_targets_{box_id}.txt"
        )
    params:
        box_id = "{box_id}"
    script:
            os.path.join("..", "..", "scripts", "process", "process_power_6_simplifier.py"
        )

rule process_simplify_all:
    input:
        [os.path.join(
        "data",
        "processed",
        "all_boxes",
        f"{box_id}",
        f"simple_network_{box_id}.gpkg") for box_id in all_boxes],
        [os.path.join(
        "data",
        "processed",
        "all_boxes",
        f"{box_id}",
        f"collapsed_sources_targets_{box_id}.txt") for box_id in all_boxes]
