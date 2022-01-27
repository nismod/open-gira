"""Process gridfinder elements for each box

"""

out_network = expand([os.path.join("data","processed","all_boxes", "{box_id}", "plants_{box_id}.gpkg"),
        os.path.join("data","processed", "all_boxes", "{box_id}", "network_{box_id}.gpkg"),
        os.path.join("data","processed", "all_boxes", "{box_id}", "targets_{box_id}.gpkg")], box_id=all_boxes),


rule process_network:
    input:
        os.path.join("data","processed", "all_boxes", "{box_id}", "powerplants_{box_id}.csv"),
        os.path.join("data","processed", "all_boxes", "{box_id}", "targets_{box_id}.csv"),
        os.path.join('data', 'processed', 'all_boxes', "{box_id}", 'gridfinder_{box_id}.gpkg'),
        os.path.join('data', 'processed', 'world_boxes_metadata.txt'),
    output:
        os.path.join("data","processed","all_boxes", "{box_id}", "plants_{box_id}.gpkg"),
        os.path.join("data","processed", "all_boxes", "{box_id}", "network_{box_id}.gpkg"),
        os.path.join("data","processed", "all_boxes", "{box_id}", "targets_{box_id}.gpkg"),
    shell:
        "python3 "+os.path.join("workflow", 'scripts', 'processing', 'process_power_4_network.py')+" {wildcards.box_id}"
