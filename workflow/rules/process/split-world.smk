"""Split the world into boxes

"""
all_box_geoms = expand(
    os.path.join(DATA_DIR, "processed", "all_boxes", "{box_id}", "geom_{box_id}.gpkg"),
    box_id=all_boxes,
)


rule world_splitter:
    input:
        os.path.join(DATA_DIR, "adminboundaries", "gadm36_levels.gpkg"),
    output:
        #os.path.join(DATA_DIR, 'processed', 'world_boxes.gpkg'),
        all_box_geoms,
        os.path.join(DATA_DIR, "processed", "world_boxes_metadata.txt"),
    shell:
        (
            "python3 "
            + os.path.join("workflow", "scripts", "process", "world_split.py")
            + f" {boxlen}"
        )
