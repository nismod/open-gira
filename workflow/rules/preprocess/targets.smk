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
    conda: "../../../environment.yml"
    input:  # note also require countries which intersect each box
        # TODO: modify script to use these inputs -- atm, paths are hardcoded in script
        rules.world_splitter.output.global_boxes,
        rules.download_worldpop_all.input.population_raster_by_country,
        os.path.join(config['output_dir'], "power_processed", "all_boxes", "{box_id}", "geom_{box_id}.gpkg"),
        os.path.join(config['output_dir'], "input", "gridfinder", "targets.tif"),
        rules.world_splitter.output.global_metadata,
        os.path.join(config['output_dir'], "input", "GDP", "GDP_per_capita_PPP_1990_2015_v2.nc"),
        os.path.join(config['output_dir'], "power_processed", "exclude_countries.json"),
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
