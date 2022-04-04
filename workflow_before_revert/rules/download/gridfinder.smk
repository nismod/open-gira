"""Download Gridfinder electricity transmission network

Reference
---------
https://gridfinder.org/
"""

out_gridfinder = [
        os.path.join('data', "gridfinder", "grid.gpkg"),
        os.path.join('data', "gridfinder", "targets.tif"),
]

# TODO follow config data_dir in script
rule download_gridfinder:
    output:
        out_gridfinder
    shell:
        """
        mkdir -p data/gridfinder
        cd data/gridfinder
        zenodo_get 10.5281/zenodo.3628142
        """
