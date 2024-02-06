"""
Intersect a network representation with hazard rasters
"""

rule rasterise_electricity_grid:
    """
    Split electricity network edges on raster grid
    Assign raster indicies to edges
    """
    input:
        network="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/edges.geoparquet",
        tif_paths=["{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/storms/wind_grid.tiff"],
    params:
        copy_raster_values=False,
    output:
        geoparquet="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/exposure/edges_split.geoparquet",
    script:
        "../../scripts/intersection.py"

"""
Test with:
snakemake --cores 1 results/power/by_country/PRI/exposure/edges_split.geoparquet
"""
