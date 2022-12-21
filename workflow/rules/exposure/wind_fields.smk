"""
Estimate max wind speed at infrastructure asset locations per event
"""


rule create_wind_grid:
    """
    Create an empty TIFF file for a given box specifying the spatial grid to
    evaluate wind speed on
    """
    input:
        box_grid="{OUTPUT_DIR}/power/world_boxes.geoparquet",
    output:
        wind_grid="{OUTPUT_DIR}/power/slice/{BOX}/storms/wind_grid.tiff",
    run:
        import os
        import geopandas as gpd

        boxes = gpd.read_parquet(input.box_grid).set_index("box_id")
        box = boxes.loc[f"box_{wildcards.BOX}", "geometry"]
        minx, miny, maxx, maxy = box.bounds
        l = config["wind_deg"]  # side length of wind cell in decimal degrees
        i = (maxx - minx) / l  # pixels in x
        j = (maxy - miny) / l  # pixels in y

        os.makedirs(os.path.dirname(output.wind_grid), exist_ok=True)
        command = f"gdal_create -outsize {i} {j} -a_srs EPSG:4326 -a_ullr {minx} {miny} {maxx} {maxy} {output.wind_grid}"
        os.system(command)

"""
Test with:
snakemake --cores 1 results/power/slice/1030/storms/wind_grid.tiff
"""


rule estimate_wind_fields:
    """
    Find maximum windspeeds for each storm for each grid cell.

    Optionally plot wind fields and save to disk
    """
    conda: "../../../environment.yml"
    input:
        storm_file = "{OUTPUT_DIR}/power/slice/{BOX}/storms/{STORM_DATASET}/tracks.geoparquet",  # TODO limited to specific storms if any
        wind_grid = "{OUTPUT_DIR}/power/slice/{BOX}/storms/wind_grid.tiff",
    output:
        # can disable plotting by setting `plot_wind_fields` to false in config
        plot_dir = directory("{OUTPUT_DIR}/power/slice/{BOX}/storms/{STORM_DATASET}/plots/"),
        wind_speeds = "{OUTPUT_DIR}/power/slice/{BOX}/storms/{STORM_DATASET}/max_wind_field.nc",
    script:
        "../../scripts/intersect/estimate_wind_fields.py"

"""
To test:
snakemake -c1 results/power/slice/1030/storms/IBTrACS/max_wind_field.nc
"""
