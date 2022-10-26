"""
Split the world into boxes
"""


rule world_splitter:
    conda: "../../../environment.yml"
    input:
        admin_data = rules.download_gadm_levels.output.admin_bounds_global_layer_per_level,
    output:
        global_metadata = "{OUTPUT_DIR}/power_processed/world_boxes_metadata.json",
        global_boxes = "{OUTPUT_DIR}/power_processed/world_boxes.gpkg",
    script:
        "../../scripts/process/world_split.py"
