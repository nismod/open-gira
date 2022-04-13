"""Finds the units that contain infrastructure

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
    params:
        region="{region}",
    script:
        os.path.join("..", "..", "scripts", "intersect", "intersect_2_gridmaker.py")


rule intersect_grid:
    input:
        region_grid,
