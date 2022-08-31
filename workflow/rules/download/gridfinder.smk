"""Download Gridfinder electricity transmission network

Reference
---------
https://gridfinder.org/
"""

out_gridfinder_global = os.path.join(config["output_dir"], "input", "gridfinder", "grid.gpkg")


rule download_gridfinder:
    output:
        out_gridfinder_global,
    shell:
        f"""
        mkdir -p {config['output_dir']}/input/gridfinder
        cd {config['output_dir']}/input/gridfinder
        zenodo_get 10.5281/zenodo.3628142
        """
