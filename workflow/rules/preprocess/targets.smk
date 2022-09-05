"""
Process targets for each country box
"""


rule process_targets:
    input:
        expand(
            os.path.join(
                config["output_dir"],
                "power_processed",
                "all_boxes",
                "{box_id}",
                "targets_{box_id}.csv",
            ),
            box_id=ALL_BOXES,
        )


rule process_target_box:
    input:  # note also require countries which intersect each box
        os.path.join(config['output_dir'], "power_processed", 'world_boxes.gpkg'),
        rules.download_population_all.input.population_raster_by_country,
        os.path.join(config['output_dir'], "power_processed", "all_boxes", "{box_id}", "geom_{box_id}.gpkg"),
        os.path.join(config['output_dir'], "input", "gridfinder", "targets.tif"),
        os.path.join(config['output_dir'], "power_processed", "world_boxes_metadata.txt"),
        os.path.join(config['output_dir'], "input", "GDP", "GDP_per_capita_PPP_1990_2015_v2.nc"),
        os.path.join(config['output_dir'], "power_processed", "exclude_countries.txt"),
    output:
        os.path.join(
            config["output_dir"],
            "power_processed",
            "all_boxes",
            "{box_id}",
            "targets_{box_id}.csv",
        ),
    params:
        box_id="{box_id}",
        output_dir=config["output_dir"],
    script:
        os.path.join("..", "..", "scripts", "process", "process_power_2_targets.py")
