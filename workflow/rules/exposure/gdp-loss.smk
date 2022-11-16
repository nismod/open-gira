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
    script:
        "../../scripts/intersect/intersect_4_gdploss.py"
