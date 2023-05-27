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
        targets = polygonise_targets(input.targets, bbox)
        targets.to_parquet(output.targets)

"""
Test with:
snakemake -c1 results/power/target_polygons.geoparquet
"""


rule annotate_targets:
    """
    Annotate targets (electricity consuming areas) with population and GDP data
    """
    conda: "../../../environment.yml"
    input:
        admin="{OUTPUT_DIR}/input/admin-boundaries/admin-level-0.geoparquet",
        population="{OUTPUT_DIR}/input/ghsl/GHS_POP_E2020_GLOBE_R2022A_54009_1000_V1_0.tif",
        targets=rules.polygonise_targets.output.targets,
        gdp="{OUTPUT_DIR}/input/GDP/GDP_per_capita_PPP_1990_2015_v2.nc",
    threads:
        config["processes_per_parallel_job"]
    output:
        targets="{OUTPUT_DIR}/power/targets.geoparquet",
    script:
        "../../scripts/preprocess/annotate_targets.py"

"""
Test with:
snakemake -c1 results/power/targets.geoparquet
"""
