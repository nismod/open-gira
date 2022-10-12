"""
Calculates the gdp losses through the storm damage.
"""

rule intersect_damages:
    conda: "../../../environment.yml"
    input:
        os.path.join(
            config["output_dir"],
            "power_intersection",
            "storm_data",
            "all_winds",
            "{storm_gcm}",
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
        CONNECTOR_OUT,
        os.path.join(
            config["output_dir"],
            "power_intersection",
            "storm_data",
            "all_winds",
            "{region}",
            "{sample}",
            "completed.txt",
        ),
    output:
        [
            os.path.join(
                config["output_dir"],
                "power_intersection",
                "storm_data",
                "individual_storms",
                "{storm_gcm}",
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
        storm_file=get_storm_file,
        all_boxes=ALL_BOXES,
    script:
        os.path.join("..", "..", "scripts", "intersect", "intersect_4_gdploss.py")
