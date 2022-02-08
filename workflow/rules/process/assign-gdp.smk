"""Assign elements for each box

"""

out_assigngdp = (
    expand(
        [
            os.path.join(
                "data",
                "processed",
                "all_boxes",
                "{box_id}",
                "network_with_gdp_{box_id}.gpkg",
            ),
            os.path.join(
                "data",
                "processed",
                "all_boxes",
                "{box_id}",
                "target_source_allocation_{box_id}.csv",
            ),
            os.path.join(
                "data",
                "processed",
                "all_boxes",
                "{box_id}",
                "edge_gdp_sorted_{box_id}.txt",
            ),
            os.path.join(
                "data",
                "processed",
                "all_boxes",
                "{box_id}",
                "targets_with_allocation_{box_id}.gpkg",
            ),
        ],
        box_id=all_boxes,
    ),
)


rule process_assigngdp:
    input:
        os.path.join(
            "data", "processed", "all_boxes", "{box_id}", "network_{box_id}.gpkg"
        ),
        os.path.join(
            "data", "processed", "all_boxes", "{box_id}", "connector_{box_id}.txt"
        ),
    output:
        os.path.join(
            "data",
            "processed",
            "all_boxes",
            "{box_id}",
            "network_with_gdp_{box_id}.gpkg",
        ),
        os.path.join(
            "data",
            "processed",
            "all_boxes",
            "{box_id}",
            "target_source_allocation_{box_id}.csv",
        ),
        os.path.join(
            "data",
            "processed",
            "all_boxes",
            "{box_id}",
            "edge_gdp_sorted_{box_id}.txt",
        ),
        os.path.join(
            "data",
            "processed",
            "all_boxes",
            "{box_id}",
            "targets_with_allocation_{box_id}.gpkg",
        ),
    shell:
        (
            "python3 "
            + os.path.join(
                "workflow", "scripts", "process", "process_power_6_assigngdp.py"
            )
            + " {wildcards.box_id}"
        )


rule process_assigngdp_all:
    input:
        out_assigngdp,
