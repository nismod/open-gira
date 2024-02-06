"""
Intersect a network representation with hazard rasters
"""


rule rasterise_osm_network:
    """
    Intersect an OSM derived network with a stack of raster maps
    """
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
