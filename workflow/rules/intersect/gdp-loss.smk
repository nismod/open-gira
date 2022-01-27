"""Calculates the gdp losses through the storm damage.

"""
import os


rule intersection_gdploss:
    input:
        os.path.join(
            "data",
            "intersection",
            "storm_data",
            "all_winds",
            "log",
            "__winds_completed_r{region}_s{sample}_y{year}.txt",  # TODO issue -> will redo all years
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
                f"edge_gdp_sorted_{box_id}.txt",
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
            "damages",
            "storm_r{region}_s{sample}_y{year}.txt",
        ),
    shell:
        (
            "python3 "
            + os.path.join("workflow", "scripts", "intersect", "intersect_4_gdploss.py")
            + " {wildcards.region} {wildcards.sample} {wildcards.year} {operationfind}"
        )
