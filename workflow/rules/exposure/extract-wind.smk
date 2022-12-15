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
        # can disable plotting by setting `plot_wind_fields` to false in config
        plot_dir = directory("{OUTPUT_DIR}/power/slice/{BOX}/storms/{STORM_DATASET}/plots/"),
        wind_speeds = "{OUTPUT_DIR}/power/slice/{BOX}/exposure/{STORM_DATASET}.nc",
    script:
        "../../scripts/intersect/intersect_3_wind_extracter.py"
