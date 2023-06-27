"""
Intersect a network representation with hazard rasters
"""


rule rasterise_osm_network:
    """
    Intersect an OSM derived network with a stack of raster maps
    """
    conda: "../../../environment.yml"
    input:
        network="{OUTPUT_DIR}/geoparquet/{DATASET}_{FILTER_SLUG}/processed/{SLICE_SLUG}_edges.geoparquet",
        tif_paths=rules.trim_hazard_data.input.trimmed_rasters,
    params:
        copy_raster_values=True
    output:
        geoparquet="{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/{SLICE_SLUG}.geoparquet",
        parquet="{OUTPUT_DIR}/splits/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/{SLICE_SLUG}.parquet",
    # intersection.py is intermittently failing, use retries as temporary fix
    # TODO: investigate why this job is sometimes failing with a coredump from intersection.py
    retries: 3
    script:
        "../../scripts/intersection.py"

"""
Test with:
snakemake --cores all results/splits/tanzania-mini_filter-road/hazard-aqueduct-river/slice-0.geoparquet
"""


rule rasterise_electricity_grid:
    """
    Split electricity network edges on raster grid
    Assign raster indicies to edges
    """
    conda: "../../../environment.yml"
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
