rule intersect_region_with_boxes:
    """
    Determine the boxes within a given region
    """
    input:
        global_boxes_metadata = "{OUTPUT_DIR}/power_processed/world_boxes_metadata.json",
        global_boxes = "{OUTPUT_DIR}/power_processed/world_boxes.gpkg",
        basin_geometry = rules.download_storm_basin_geometry.output.geometry
    params:
        region_name = "{wildcards.STORM_BASIN}"
    output:
        boxes_in_region = "{OUTPUT_DIR}/power_intersection/regions/{STORM_BASIN}_boxes.json",
        boxes_in_region_with_assets = "{OUTPUT_DIR}/power_intersection/regions/{STORM_BASIN}_boxes_with_assets.json"
    script:
        os.path.join("..", "..", "scripts", "intersect", "intersect_1_regions.py")
