"""Download Gridfinder electricity transmission network

Reference
---------
https://gridfinder.org/
"""


ELECTRICITY_GRID_GLOBAL = os.path.join(config["output_dir"], "input", "gridfinder", "grid.gpkg")


rule download_gridfinder:
    output:
        ELECTRICITY_GRID_GLOBAL,
    shell:
        f"""
        mkdir -p {config['output_dir']}/input/gridfinder
        cd {config['output_dir']}/input/gridfinder
        zenodo_get 10.5281/zenodo.3628142
        """
