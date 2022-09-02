"""
Finds the units that contain infrastructure
"""

import os


rule intersect_unit_generator:
    input:
        expand(
            os.path.join(
                config["output_dir"],
                "power_processed",
                "all_boxes",
                "{box_id}",
                "geom_{box_id}.gpkg",
            ),
            box_id=ALL_BOXES,
        ),
        os.path.join(
            config["output_dir"],
            "input",
            "stormtracks",
            "fixed",
            "STORM_FIXED_RETURN_PERIODS_{region}.nc",
        ),
        os.path.join(
            config["output_dir"], "power_intersection", "regions", "{region}_boxes_with_assets.json"
        ),
    output:
        os.path.join(
            config["output_dir"], "power_intersection", "regions", "{region}_unit.gpkg"
        ),
    params:
        region="{region}",
        output_dir=config["output_dir"],
    script:
        os.path.join("..", "..", "scripts", "intersect", "intersect_2_gridmaker.py")


rule intersect_grid:
    input:
        expand(
            os.path.join(
                config["output_dir"], "power_intersection", "regions", "{region}_unit.gpkg"
            ),
            region=REGIONS,
        )
