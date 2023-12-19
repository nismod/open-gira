rule simplify_admin_bounds:
    input:
        all_admin_bounds = rules.download_gadm_levels.output.gpkg
    output:
        simplified_layer = "{OUTPUT_DIR}/input/admin-boundaries/{ADMIN_SLUG}.geoparquet"
    run:
        import geopandas as gpd

        df = gpd.read_file(input.all_admin_bounds, layer=int(wildcards.ADMIN_SLUG.replace("admin-level-", "")))

        TOLERANCE_METRES = 50
        # GADM geometry is very precise -- drop some precision to save on memory and disk
        # https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.simplify.html
        original_crs = df.crs
        df["geometry"] = df["geometry"].to_crs(epsg=3857).simplify(TOLERANCE_METRES).to_crs(crs=original_crs)

        df.to_parquet(output.simplified_layer)

"""
To test:
snakemake --cores 1 results/input/admin-boundaries/admin-level-0.geoparquet
"""
