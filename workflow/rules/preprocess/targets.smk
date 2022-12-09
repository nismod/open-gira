"""
Process targets to estimate population and GDP
"""

rule process_target_box:
    conda: "../../../environment.yml"
    input:
        population="{OUTPUT_DIR}/input/ghsl/GHS_POP_E2020_GLOBE_R2022A_54009_1000_V1_0.tif",
        targets=rules.download_gridfinder.output.electricity_targets_global,
        gdp="{OUTPUT_DIR}/input/GDP/GDP_per_capita_PPP_1990_2015_v2.nc",
        global_boxes=rules.world_splitter.output.global_boxes,
    output:
        targets="{OUTPUT_DIR}/power/slice/{BOX}/network/targets_{BOX}.geoparquet",
    script:
        "../../scripts/process/process_power_2_targets.py"
