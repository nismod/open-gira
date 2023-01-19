"""Download Gridfinder electricity transmission network

Reference
---------
https://gridfinder.org/
"""


rule download_gridfinder:
    conda: "../../../environment.yml"
    output:
        electricity_grid_global = "{OUTPUT_DIR}/input/gridfinder/grid.gpkg",
        electricity_targets_global = "{OUTPUT_DIR}/input/gridfinder/targets.tif",
    shell:
        f"""
        mkdir -p {config['output_dir']}/input/gridfinder
        cd {config['output_dir']}/input/gridfinder
        zenodo_get 10.5281/zenodo.3628142
        """

"""
To test:
snakemake --cores 1 results/input/gridfinder/grid.gpkg
"""


rule gridfinder_to_geoparquet:
    input:
        geopackage = rules.download_gridfinder.output.electricity_grid_global,
    output:
        geoparquet = "{OUTPUT_DIR}/power/gridfinder.geoparquet",
    run:
        import geopandas as gpd
        gpd.read_file(input.geopackage).to_parquet(output.geoparquet)

"""
To test:
snakemake --cores 1 results/input/gridfinder/grid.geoparquet
"""
