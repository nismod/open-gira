"""
Process gridfinder elements for each box
"""


rule process_gridfinder:
    conda: "../../../environment.yml"
    input:
        gridfinder="{OUTPUT_DIR}/input/gridfinder/grid.geoparquet",
        global_boxes=rules.world_splitter.output.global_boxes,
    output:
        gridfinder="{OUTPUT_DIR}/processed/power/{BOX}/gridfinder_{BOX}.geoparquet",
    script:
       "../../scripts/process/process_power_3_gridfinder.py"
