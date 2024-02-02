rule polygonise_targets:
    """
    Process raster of target locations into set of polygons

    Running time: ~5min
    """
    input:
        targets=rules.download_gridfinder.output.electricity_targets_global,
    output:
        targets="{OUTPUT_DIR}/power/target_polygons.geoparquet",
    run:
        from open_gira.grid import polygonise_targets

        import shapely
        from shapely.geometry.polygon import Polygon
        from shapely.geometry import box

        # do not process polar regions
        bbox: Polygon = shapely.geometry.box(-180, -60, 180, 60)

        # manually set CRS of raster (and returned targets) using EPSG code
        # allowing polygonise_targets to infer it from the raster metadata
        # results in geopandas complaining the CRSs are not equal when their
        # transforms are in fact equal, but one is a proj string, another an
        # EPSG code
        targets = polygonise_targets(input.targets, bbox, 4326)
        targets.to_parquet(output.targets)

"""
Test with:
snakemake -c1 results/power/target_polygons.geoparquet
"""


rule reproject_targets_to_mollweide:
    """
    For calculating the population of each target, first transform to population raster CRS.
    """
    input:
        targets=rules.polygonise_targets.output.targets,
    output:
        # use gpkg so exactextract (GDAL) can read
        targets="{OUTPUT_DIR}/power/target_polygons_mollweide.gpkg",
    run:
        import geopandas as gpd

        targets = gpd.read_parquet(input.targets)
        targets.to_crs("+proj=moll").to_file(output.targets)

"""
Test with:
snakemake -c1 results/power/target_polygons_mollweide.gpkg
"""


rule calculate_population_of_targets:
    """
    Use exactextract to calculate the population sum overlapping with each target polygon.
    """
    input:
        polygons=rules.reproject_targets_to_mollweide.output.targets,
        raster="{OUTPUT_DIR}/input/ghsl/GHS_POP_E2020_GLOBE_R2022A_54009_1000_V1_0.tif",
    output:
        population_table="{OUTPUT_DIR}/power/target_population.csv",
    shell:
        """
        exactextract \
            --polygons {input.polygons} \
            --raster "population:{input.raster}" \
            --stat "sum(population)" \
            --output {output.population_table} \
            --fid "id"
        """

"""
Test with:
snakemake -c1 results/power/target_population.csv
"""


rule annotate_targets:
    """
    Annotate targets (electricity consuming areas) with population and GDP data
    """
    input:
        admin="{OUTPUT_DIR}/input/admin-boundaries/admin-level-0.geoparquet",
        population=rules.calculate_population_of_targets.output.population_table,
        targets=rules.polygonise_targets.output.targets,
        gdp=rules.download_GDP.output.gdp_pc,
    output:
        targets="{OUTPUT_DIR}/power/targets.geoparquet",
    script:
        "./annotate_targets.py"

"""
Test with:
snakemake -c1 results/power/targets.geoparquet
"""


checkpoint rank_countries_by_target_count:
    """
    Create a lookup table from ISO alpha-3 country code to number of targets in
    that state.

    This is used for setting resource limits (cores) for analysis jobs (e.g.
    wind field estimation).
    """
    input:
        targets = rules.annotate_targets.output.targets
    output:
        lookup_table = "{OUTPUT_DIR}/power/target_count_by_country.csv"
    run:
        import geopandas as gpd

        targets = gpd.read_parquet(input.targets)
        target_rank = targets.value_counts("iso_a3", ascending=True).reset_index().rename(columns={"count": "target_count"})
        target_rank.set_index("iso_a3")
        target_rank.to_csv(output.lookup_table)

"""
Test with:
snakemake -c1 results/power/target_count_by_country.csv
"""
