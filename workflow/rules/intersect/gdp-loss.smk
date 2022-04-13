"""Calculates the gdp losses through the storm damage.

"""


rule intersection_gdp_loss:
    input:
        os.path.join(
            "data",
            "intersection",
            "storm_data",
            "all_winds",
            "{region}",
            "{sample}",
            "TC_r{region}_s{sample}_n{nh}.csv",
        ),
        os.path.join("data", "intersection", "regions", "{region}_unit.gpkg"),
        [
            os.path.join(
                "data", "processed", "all_boxes", f"{box_id}", f"targets_{box_id}.gpkg"
            )
            for box_id in all_boxes
        ],
        out_connector,
        os.path.join(
            "data",
            "intersection",
            "storm_data",
            "all_winds",
            "{region}",
            "{sample}",
            "{region}_{sample}_completed.txt",
        ),
    output:
        os.path.join(
            "data",
            "intersection",
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
    script:
        os.path.join("..", "..", "scripts", "intersect", "intersect_4_gdploss.py")
