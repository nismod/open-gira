"""Download Gridfinder electricity transmission network

Reference
---------
https://gridfinder.org/
"""

# TODO follow config data_dir in script
rule download_gridfinder:
    output:
        os.path.join(config['data_dir'], "gridfinder", "grid.gpkg"),
        os.path.join(config['data_dir'], "gridfinder", "targets.tif"),,
    shell:
        """
        mkdir -p data/gridfinder
        cd data/gridfinder
        zenodo_get 10.5281/zenodo.3628142
        """
