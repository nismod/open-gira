"""
Process targets to estimate population and GDP
"""

rule process_target_box:
    conda: "../../../environment.yml"
    input:
        population="{OUTPUT_DIR}/input/ghsl/GHS_POP_E2020_GLOBE_R2022A_54009_1000_V1_0.tif",
        targets="{OUTPUT_DIR}/input/gridfinder/targets.tif",
        gdp="{OUTPUT_DIR}/input/GDP/GDP_per_capita_PPP_1990_2015_v2.nc",
    output:
        targets="{OUTPUT_DIR}/processed/power_targets.geoparquet",
    script:
        "../../scripts/process/process_power_2_targets.py"
