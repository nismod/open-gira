"""Process targets for each country box

"""


out_targets = expand(
    os.path.join("data", "processed", "all_boxes", "{box_id}", "targets_{box_id}.csv"),
    box_id=all_boxes,
)


rule process_targets:
    input:
        out_targets,


rule process_target_box:
    input:  # note also require countries which intersect each box
        #os.path.join('data', 'processed', 'world_boxes.gpkg'),
        os.path.join("data", "processed", "all_boxes", "{box_id}", "geom_{box_id}.gpkg"),
        os.path.join("data", "gridfinder", "targets.tif"),
        os.path.join("data", "processed", "world_boxes_metadata.txt"),
        os.path.join("data", "GDP", "GDP_per_capita_PPP_1990_2015_v2.nc"),
        os.path.join("data", "adminboundaries", "exclude_countries.txt"),
    output:
        os.path.join(
            "data", "processed", "all_boxes", "{box_id}", "targets_{box_id}.csv"
        ),
    params:
        box_id="{box_id}",
    script:
        os.path.join("..", "..", "scripts", "process", "process_power_2_targets.py")
