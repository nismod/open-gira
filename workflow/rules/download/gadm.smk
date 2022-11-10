"""
Download GADM boundaries

Reference
---------
https://gadm.org/data.html
"""


rule download_gadm_levels:
    output:
        admin_bounds_global_layer_per_level = os.path.join(
            config['output_dir'], "input", "admin-boundaries", "gadm36_levels.gpkg"
        )
    shell:
        f"""
        wget https://geodata.ucdavis.edu/gadm/gadm3.6/gadm36_levels_gpkg.zip \
            --output-document={config['output_dir']}/input/admin-boundaries/gadm36_levels_gpkg.zip
        unzip -o {config['output_dir']}/input/admin-boundaries/gadm36_levels_gpkg.zip \
            -d {config['output_dir']}/input/admin-boundaries
        rm {config['output_dir']}/input/admin-boundaries/gadm36_levels_gpkg.zip
        """
