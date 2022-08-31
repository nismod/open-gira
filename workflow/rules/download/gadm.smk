"""
Download GADM boundaries

Reference
---------
https://gadm.org/data.html
"""

import os


rule download_gadm:
    output:
        ADMIN_BOUNDS_GLOBAL_SINGLE_LAYER,
    shell:
        f"""
        wget https://geodata.ucdavis.edu/gadm/gadm3.6/gadm36_gpkg.zip \
            --directory-prefix={config['output_dir']}/input/admin-boundaries
        unzip -o {config['output_dir']}/input/admin-boundaries/gadm36_gpkg.zip -d {config['output_dir']}/input/admin-boundaries
        """


"""
Test with:
snakemake --cores 1 results/input/admin-boundaries/gadm36.gpkg
"""


rule download_gadm_levels:
    output:
        ADMIN_BOUNDS_GLOBAL_LAYER_PER_LEVEL,
    shell:
        f"""
        wget https://geodata.ucdavis.edu/gadm/gadm3.6/gadm36_levels_gpkg.zip \
            --output-document={config['output_dir']}/input/admin-boundaries/gadm36_levels_gpkg.zip
        unzip -o {config['output_dir']}/input/admin-boundaries/gadm36_levels_gpkg.zip -d {config['output_dir']}/input/admin-boundaries
        """


# download admin boundaries per country
rule download_gadm_by_country:
    output:
        os.path.join(config['output_dir'], "input", "admin-boundaries", "gadm36_{code}.gpkg"),
    shell:
        f"""
        wget https://geodata.ucdavis.edu/gadm/gadm3.6/gpkg/gadm36_{{wildcards.code}}_gpkg.zip \
            --output-document={config['output_dir']}/input/admin-boundaries/gadm36_{{wildcards.code}}_gpkg.zip
        unzip -o {config['output_dir']}/input/admin-boundaries/gadm36_{{wildcards.code}}_gpkg.zip -d {config['output_dir']}/input/admin-boundaries
        """
