"""
Intersect a network representation with hazard rasters
"""


rule network_intersection:
    """
    Intersect network nodes and edges with a stack of raster maps
    """
    input:
        nodes="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/nodes.geoparquet",
        edges="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/network/edges.geoparquet",
        tif_paths=rules.trim_hazard_data.input.trimmed_rasters,
    output:
        nodes="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/{HAZARD_SLUG}_{DATASET}/nodes.geoparquet",
        edges="{OUTPUT_DIR}/power/by_country/{COUNTRY_ISO_A3}/{HAZARD_SLUG}_{DATASET}/edges.geoparquet",
    run:
        import geopandas
        from snail.intersection import (
            prepare_linestrings,
            prepare_points,
            split_points,
            split_linestrings,
            split_features_for_rasters,
        )
        from snail.io import (
            associate_raster_files,
            extend_rasters_metadata,
        )
        from open_gira import fields

        raster_paths: list[str] = input.tif_paths
        raster_keys = [
            f"{fields.HAZARD_PREFIX}{Path(path).name}"
            for path in raster_paths
        ]
        breakpoint()

        rasters, transforms = extend_rasters_metadata(rasters)

        nodes = geopandas.read_parquet(input.nodes)
        nodes_prepared = prepare_points(nodes)
        nodes_split = split_features_for_rasters(nodes_prepared, transforms, split_points)
        nodes_with_data = associate_raster_files(nodes_split, rasters)
        nodes_with_data.to_parquet(output.nodes)

        edges = geopandas.read_parquet(input.edges)
        edges_prepared = prepare_linestrings(edges)
        edges_split = split_features_for_rasters(edges_prepared, transforms, split_linestrings)
        edges_with_data = associate_raster_files(edges_split, rasters)
        edges_with_data.to_parquet(output.edges)

"""
Test with:
snakemake --cores all results/power/by_country/TZA/hazard-aqueduct-river_planet-latest/nodes.geoparquet
"""
