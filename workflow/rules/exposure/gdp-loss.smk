"""
Calculates the gdp losses through the storm damage.
"""

rule intersect_damages:
    conda: "../../../environment.yml"
    input:
        # do we need all boxes assessed? Yes, think so
        edges = expand("{{OUTPUT_DIR}}/processed/power/{box}/edges_{box}.parquet", box=ALL_BOXES),
        nodes = expand("{{OUTPUT_DIR}}/processed/power/{box}/nodes_{box}.parquet", box=ALL_BOXES),
        wind_speeds = expand(
            "{{OUTPUT_DIR}}/processed/power/exposure/{STORM_SUBSET}/{STORM_MODEL}/{BASIN}/ws_{BOX}_{SAMPLE}.parquet",
            box=ALL_BOXES
        ),
    output:
        # rows are "targets" - areas served by the power network
        # columns are target_id, storm_id, fraction_of_service, customers_affected (population * f), gdp_loss (gdp_pc * customers_affected)
        losses = "{OUTPUT_DIR}/processed/power/damages/{STORM_SUBSET}/{STORM_MODEL}/{BASIN}/ws_{SAMPLE}_{THRESHOLD}.parquet"
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
        all_boxes=ALL_BOXES,
    script:
        "../../scripts/intersect/intersect_4_gdploss.py"
