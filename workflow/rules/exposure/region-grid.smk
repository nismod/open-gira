"""
Intersect power network with a wind speed grid

- Split network edges on raster grid
- Assign raster indicies to edges
"""


rule create_wind_grid:
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


rule rasterise_network:
    conda: "../../../environment.yml"
    input:
        network="{OUTPUT_DIR}/power/slice/{BOX}/network/edges_{BOX}.geoparquet",
        tif_paths=["{OUTPUT_DIR}/power/slice/{BOX}/storms/wind_grid.tiff"],
    params:
        copy_raster_values=False,
    output:
        geoparquet="{OUTPUT_DIR}/power/slice/{BOX}/exposure/edges_split_{BOX}.geoparquet",
    script:
        "../../scripts/intersection.py"

"""
Test with:
snakemake --cores 1 results/power/slice/1030/exposure/edges_split_1030.geoparquet
"""
