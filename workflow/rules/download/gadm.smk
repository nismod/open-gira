"""
Download GADM boundaries

Reference
---------
https://gadm.org/data.html
"""


rule download_gadm:
    """
    Test with:
    snakemake --cores 1 results/input/admin-boundaries/gadm36.gpkg
    """
    output:
        admin_bounds_global_single_layer = os.path.join(config['output_dir'], "input", "admin-boundaries", "gadm36.gpkg")
    shell:
        f"""
        wget https://geodata.ucdavis.edu/gadm/gadm3.6/gadm36_gpkg.zip \
            --directory-prefix={config['output_dir']}/input/admin-boundaries
        unzip -o {config['output_dir']}/input/admin-boundaries/gadm36_gpkg.zip \
            -d {config['output_dir']}/input/admin-boundaries
        rm {config['output_dir']}/input/admin-boundaries/gadm36_gpkg.zip
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


rule download_gadm_by_country:
    output:
        "{OUTPUT_DIR}/input/admin-boundaries/gadm36_{CODE}.gpkg"
    shell:
        f"""
        wget https://geodata.ucdavis.edu/gadm/gadm3.6/gpkg/gadm36_{{wildcards.CODE}}_gpkg.zip \
            --output-document={{wildcards.OUTPUT_DIR}}/input/admin-boundaries/gadm36_{{wildcards.CODE}}_gpkg.zip
        unzip -o {{wildcards.OUTPUT_DIR}}/input/admin-boundaries/gadm36_{{wildcards.CODE}}_gpkg.zip \
            -d {{wildcards.OUTPUT_DIR}}/input/admin-boundaries
        rm {{wildcards.OUTPUT_DIR}}/input/admin-boundaries/gadm36_{{wildcards.CODE}}_gpkg.zip
        """


rule download_gadm_all_countries:
    input:
        admin_bounds_file_per_country = expand(
            os.path.join(
                config["output_dir"], "input", "admin-boundaries", "gadm36_{code}.gpkg"
            ),
            code=COUNTRY_CODES,
        )
