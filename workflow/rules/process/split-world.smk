"""Split the world into boxes

"""


rule world_splitter:
    input:
        os.path.join(DATA_DIR,"adminboundaries","gadm36_levels.gpkg")
    output:
        os.path.join(DATA_DIR, 'processed', 'world_boxes.gpkg'),
        os.path.join(DATA_DIR, 'processed', 'world_boxes_metadata.txt'),
    shell:
        "python3 "+os.path.join(WORKFLOW_DIR, 'scripts', 'processing', 'world_split.py')+f" {boxlen}"