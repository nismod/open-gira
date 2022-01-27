"""Extracts the box_ids within the selected region

"""
import os

region_box = expand(
    os.path.join("data", "intersection", "regions", "{region}_boxes.txt"),
    region=REGIONS,
)


rule intersect_regions_indiv:
    input:
        os.path.join("data", "processed", "world_boxes_metadata.txt"),
        #os.path.join("data", "processed", "world_boxes.gpkg"),  # removed because opening on QGIS tampers with metadata
        os.path.join(
            "data", "stormtracks", "fixed", "STORM_FIXED_RETURN_PERIODS_{region}.nc"
        ),
        expand(
            os.path.join(
                "data",
                "processed",
                "all_boxes",
                "{box_id}",
                "gridfinder_{box_id}.gpkg",
            ),
            box_id=all_boxes,
        ),
    output:
        os.path.join("data", "intersection", "regions", "{region}_boxes.txt"),
    shell:
        (
            "python3 "
            + os.path.join("workflow", "scripts", "intersect", "intersect_1_regions.py")
            + " {wildcards.region}"
        )


rule intersect_regions:
    input:
        region_box,
