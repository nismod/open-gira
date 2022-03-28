"""Calculates the gdp losses through the storm damage.

"""


rule intersection_gdploss:
    input:
        os.path.join(
            "data",
            "intersection",
            "storm_data",
            "all_winds",
            "{region}",
            "TC_r{region}_s{sample}_n{nh}.csv",
        ),
        os.path.join("data", "intersection", "regions", "{region}_unit.gpkg"),
        [
            os.path.join(
                "data",
                "processed",
                "all_boxes",
                f"{box_id}",
                f"network_with_gdp_{box_id}.gpkg",
            )
            for box_id in all_boxes
        ],
        [
            os.path.join(
                "data",
                "processed",
                "all_boxes",
                f"{box_id}",
                f"edge_gdp_sorted_{box_id}.hkl",
            )
            for box_id in all_boxes
        ],
        [
            os.path.join(
                "data",
                "processed",
                "all_boxes",
                f"{box_id}",
                f"component_indices_{box_id}.hkl",
            )
            for box_id in all_boxes
        ],
        [
            os.path.join(
                "data", "processed", "all_boxes", f"{box_id}", f"targets_{box_id}.gpkg"
            )
            for box_id in all_boxes
        ],
    output:
        os.path.join(
            "data",
            "intersection",
            "storm_data",
            "individual_storms",
            "{region}",
            "storm_{nh}",
            "storm_r{region}_s{sample}_n{nh}.txt",
        ),
        os.path.join(
            "data",
            "intersection",
            "storm_data",
            "individual_storms",
            "{region}",
            "storm_{nh}",
            "storm_track_r{region}_s{sample}_n{nh}.gpkg",
        ),
    params:
        region = "{region}", sample = "{sample}", nh = "{nh}", op_find = operationfind,
    script:
            os.path.join("..", "..", "scripts", "intersect", "intersect_4_gdploss.py"
        )