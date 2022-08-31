"""Calculates the gdp losses through the storm damage.

"""

rule intersect_damages:
    input:
        os.path.join(
            config["output_dir"],
            "power_intersection",
            "storm_data",
            "all_winds",
            "{region}",
            "{sample}",
            "TC_r{region}_s{sample}_n{nh}.csv",
        ),
        os.path.join(
            config["output_dir"], "power_intersection", "regions", "{region}_unit.gpkg"
        ),
        [
            os.path.join(
                config["output_dir"],
                "power_processed",
                "all_boxes",
                f"{box_id}",
                f"targets_{box_id}.gpkg",
            )
            for box_id in ALL_BOXES
        ],
        out_connector,
        os.path.join(
            config["output_dir"],
            "power_intersection",
            "storm_data",
            "all_winds",
            "{region}",
            "{sample}",
            "{region}_{sample}_completed.txt",
        ),
    output:
        [
            os.path.join(
                config["output_dir"],
                "power_intersection",
                "storm_data",
                "individual_storms",
                "{region}",
                "{sample}",
                "storm_{nh}",
                f"{thrval}",
                "storm_r{region}_s{sample}_n{nh}.txt",
            )
            for thrval in WIND_SPEED_THRESHOLDS_MS
        ],
    params:
        region="{region}",
        sample="{sample}",
        nh="{nh}",
        output_dir=config["output_dir"],
        reconstruction_cost_lowmedium=config["reconstruction_cost_lowmedium"],
        reconstruction_cost_high=config["reconstruction_cost_high"],
        central_threshold=config["central_threshold"],
        minimum_threshold=config["minimum_threshold"],
        maximum_threshold=config["maximum_threshold"],
        wind_file_start=wind_file_start,
        wind_file_end=wind_file_end,
        all_boxes=ALL_BOXES,
    script:
        os.path.join("..", "..", "scripts", "intersect", "intersect_4_gdploss.py")
