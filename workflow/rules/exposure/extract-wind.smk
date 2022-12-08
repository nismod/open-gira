"""
Estimate max wind speed at infrastructure asset locations per event
"""

rule estimate_wind_fields:
    """
    Find maximum windspeeds for each storm for each grid cell.
    """
    conda: "../../../environment.yml"
    input:
        edges_split = "{OUTPUT_DIR}/power/slice/{BOX}/exposure/edges_split_{BOX}.geoparquet",
        storm_file = "{OUTPUT_DIR}/power/slice/{BOX}/storms/{STORM_DATASET}/tracks.geoparquet",  # TODO limited to specific storms if any
        wind_grid = "{OUTPUT_DIR}/power/slice/{BOX}/storms/wind_grid.tiff",
    output:
        plot_dir = "{OUTPUT_DIR}/power/slice/{BOX}/storms/{STORM_DATASET}/plots/",  # set to empty string to disable plotting
        wind_speeds = "{OUTPUT_DIR}/power/slice/{BOX}/exposure/{STORM_DATASET}.parquet",
    script:
        "../../scripts/intersect/intersect_3_wind_extracter.py"
