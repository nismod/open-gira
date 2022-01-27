"""Finds which units contain infrastructure

"""
import os

region_grid = expand(
    os.path.join("data", "intersection", "regions", "{region}_unit.gpkg"),
    region=REGIONS,
)


rule intersect_grid_indiv:
    input:
        expand(
            os.path.join(
                "data", "processed", "all_boxes", "{box_id}", "geom_{box_id}.gpkg"
            ),
            box_id=all_boxes,
        ),
        os.path.join(
            "data", "stormtracks", "fixed", "STORM_FIXED_RETURN_PERIODS_{region}.nc"
        ),
        os.path.join("data", "intersection", "regions", "{region}_boxes.txt"),
    output:
        os.path.join("data", "intersection", "regions", "{region}_unit.gpkg"),
    shell:
        (
            "python3 "
            + os.path.join(
                "workflow", "scripts", "intersect", "intersect_2_gridmaker.py"
            )
            + " {wildcards.region}"
        )


rule intersect_grid:
    input:
        region_grid,
