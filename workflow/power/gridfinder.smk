"""Download Gridfinder electricity transmission network

Reference
---------
https://gridfinder.org/
"""


rule download_gridfinder:
    output:
        electricity_grid_global = "{OUTPUT_DIR}/input/gridfinder/grid.gpkg",
        electricity_targets_global = "{OUTPUT_DIR}/input/gridfinder/targets.tif",
    shell:
        """
        mkdir -p {wildcards.OUTPUT_DIR}/input/gridfinder
        cd {wildcards.OUTPUT_DIR}/input/gridfinder
        zenodo_get 10.5281/zenodo.3628142
        """

"""
To test:
snakemake --cores 1 results/input/gridfinder/grid.gpkg
"""
