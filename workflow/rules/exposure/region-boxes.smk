rule intersect_region_with_boxes:
    """
    Determine the boxes within a given region
    """
    conda: "../../../environment.yml"
    input:
        global_boxes_metadata = rules.world_splitter.output.global_metadata,
        global_boxes = rules.world_splitter.output.global_boxes,
        basin_geometry = "bundled_data/storm_basins.geojson",
    output:
        boxes_in_region = "{OUTPUT_DIR}/processed/power/{REGION}_boxes.json",
        boxes_in_region_with_assets = "{OUTPUT_DIR}/processed/power/{REGION}_boxes_with_assets.json"
    script:
       "../../scripts/intersect/intersect_1_regions.py"
