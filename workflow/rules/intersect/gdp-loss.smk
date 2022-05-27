"""Calculates the gdp losses through the storm damage.

"""


rule intersect_damages:
    input:
        os.path.join(
            config['output_dir'],
            "power_intersection",
            "storm_data",
            "all_winds",
            "{region}",
            "{sample}",
            "TC_r{region}_s{sample}_n{nh}.csv",
        ),
        os.path.join(config['output_dir'], "power_intersection", "regions", "{region}_unit.gpkg"),
        [
            os.path.join(
                config['output_dir'], "power_processed", "all_boxes", f"{box_id}", f"targets_{box_id}.gpkg"
            )
            for box_id in all_boxes
        ],
        out_connector,
        os.path.join(
            config['output_dir'],
            "power_intersection",
            "storm_data",
            "all_winds",
            "{region}",
            "{sample}",
            "{region}_{sample}_completed.txt",
        ),
    output:
        os.path.join(
            config['output_dir'],
            "power_intersection",
            "storm_data",
            "individual_storms",
            "{region}",
            "{sample}",
            "storm_{nh}",
            "storm_r{region}_s{sample}_n{nh}.txt",
        ),
    params:
        region="{region}",
        sample="{sample}",
        nh="{nh}",
        output_dir = config['output_dir'],
        reconstruction_cost = config['reconstruction_cost']
    script:
        os.path.join("..", "..", "scripts", "intersect", "intersect_4_gdploss.py")
