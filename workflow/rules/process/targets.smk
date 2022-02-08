"""Process targets for each country box

"""

exts = ["csv"]
out_targets = expand(
    os.path.join(
        DATA_DIR, "processed", "all_boxes", "{box_id}", "targets_{box_id}.{ext}"
    ),
    box_id=all_boxes,
    ext=exts,
)


rule process_targets:
    input:
        out_targets,


rule process_target_box:
    input:  # note also require countries which intersect each box
        #os.path.join(DATA_DIR, 'processed', 'world_boxes.gpkg'),
        os.path.join(
            DATA_DIR, "processed", "all_boxes", "{box_id}", "geom_{box_id}.gpkg"
        ),
        os.path.join(DATA_DIR, "gridfinder", "targets.tif"),
        os.path.join(DATA_DIR, "processed", "world_boxes_metadata.txt"),
        os.path.join(DATA_DIR, "GDP", "GDP_per_capita_PPP_1990_2015_v2.nc"),
        os.path.join("data", "adminboundaries", "exclude_countries.txt"),
    output:
        os.path.join(
            DATA_DIR, "processed", "all_boxes", "{box_id}", "targets_{box_id}.csv"
        ),
    shell:
        (
            "python3 "
            + os.path.join(
                "workflow", "scripts", "process", "process_power_2_targets.py"
            )
            + " {wildcards.box_id}"
        )
