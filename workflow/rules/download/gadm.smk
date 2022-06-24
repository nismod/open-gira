"""Download GADM boundaries

Reference
---------
https://gadm.org/data.html
"""

out_adminboundaries = os.path.join(config['output_dir'], "input", "adminboundaries", "gadm36.gpkg")
out_adminboundaries_codes = expand(
    os.path.join(config['output_dir'], "input", "adminboundaries", "gadm36_{code}.gpkg"), code=COUNTRY_CODES
)


rule download_gadm:
    output:
        out_adminboundaries,
    shell:
        f"""
        wget https://geodata.ucdavis.edu/gadm/gadm3.6/gpkg/gadm36_gpkg.zip \
            --directory-prefix={config['output_dir']}/input/adminboundaries
        unzip -o {config['output_dir']}/input/adminboundaries/gadm36_gpkg.zip -d {config['output_dir']}/input/adminboundaries
        """


# download admin boundaries per country
rule download_gadm_by_country:
    output:
        os.path.join(config['output_dir'], "input", "adminboundaries", "gadm36_{code}.gpkg"),
    shell:
        """
        wget https://geodata.ucdavis.edu/gadm/gadm3.6/gpkg/gadm36_{wildcards.code}_gpkg.zip \
            --directory-prefix={config['output_dir']}/input/adminboundaries
        unzip -o {config['output_dir']}/input/adminboundaries/gadm36_{wildcards.code}_gpkg.zip -d {config['output_dir']}/input/adminboundaries
        """


out_adminboundaries_levels = os.path.join(
    config['output_dir'], "input", "adminboundaries", "gadm36_levels.gpkg"
)


rule download_gadm_levels:
    output:
        out_adminboundaries_levels,
    shell:
        f"""
        wget https://geodata.ucdavis.edu/gadm/gadm3.6/gadm36_levels_gpkg.zip \
            --directory-prefix={config['output_dir']}/input/adminboundaries
        unzip -o {config['output_dir']}/input/adminboundaries/gadm36_levels_gpkg.zip -d {config['output_dir']}/input/adminboundaries
        """
