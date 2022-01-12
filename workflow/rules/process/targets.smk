"""Process targets for each country box

"""

exts = ['csv']
out_targets = expand(os.path.join(DATA_DIR, 'processed', 'all_boxes', "{box_id}", "targets_{box_id}.{ext}"), box_id=all_boxes, ext=exts)


#out_targets = [os.path.join(DATA_DIR, 'processed', 'all_boxes', f"box_{idx}", f"targets_box_{idx}.csv") for idx in [1940, 1941,1942]]  # TODO temp

rule process_targets:
    input:
        out_targets


rule process_target_box:
    input:  # note also require countries which intersect each box
        os.path.join(DATA_DIR, 'processed', 'world_boxes.gpkg'),
        os.path.join(DATA_DIR,"gridfinder","targets.tif"),
        os.path.join(DATA_DIR, 'processed', 'world_boxes_metadata.txt'),
        os.path.join(DATA_DIR,"GDP","GDP_per_capita_PPP_1990_2015_v2.nc"),
    output:
        os.path.join(DATA_DIR, 'processed', 'all_boxes', "{box_id}", "targets_{box_id}.csv"),
    shell:
        "python3 "+os.path.join(WORKFLOW_DIR, 'scripts', 'processing', 'process_power_2_targets.py')+" {wildcards.box_id}"