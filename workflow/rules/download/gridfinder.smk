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
