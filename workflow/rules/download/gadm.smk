"""
Download GADM boundaries

Reference
---------
https://gadm.org/data.html
"""


rule download_gadm_levels:
    output:
        gpkg = "{OUTPUT_DIR}/input/admin-boundaries/gadm36_levels.gpkg"
    shell:
        """
        wget https://geodata.ucdavis.edu/gadm/gadm3.6/gadm36_levels_gpkg.zip \
            --output-document={wildcards.OUTPUT_DIR}/input/admin-boundaries/gadm36_levels_gpkg.zip
        unzip -o {wildcards.OUTPUT_DIR}/input/admin-boundaries/gadm36_levels_gpkg.zip \
            -d {wildcards.OUTPUT_DIR}/input/admin-boundaries
        rm {wildcards.OUTPUT_DIR}/input/admin-boundaries/gadm36_levels_gpkg.zip
        """

"""
Test with:
snakemake -c1 -- results/input/admin-boundaries/gadm36_levels.gpkg
"""