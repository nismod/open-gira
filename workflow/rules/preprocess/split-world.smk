"""
Split the world into boxes
"""


rule world_splitter:
    conda: "../../../environment.yml"
    input:
        admin_data = "{OUTPUT_DIR}/input/admin-boundaries/admin-level-0.geoparquet",
    output:
        global_metadata = "{OUTPUT_DIR}/power/world_boxes_metadata.json",
        global_boxes = "{OUTPUT_DIR}/power/world_boxes.geoparquet",
    script:
        "../../scripts/preprocess/world_split.py"
