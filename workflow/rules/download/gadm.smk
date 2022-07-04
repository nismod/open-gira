"""Download GADM boundaries

Reference
---------
https://gadm.org/data.html
"""
import os


out_adminboundaries_codes = expand(
    os.path.join(config['output_dir'], "input", "adminboundaries", "gadm36_{code}.gpkg"), code=COUNTRY_CODES
)


rule download_gadm:
    output:
        "{OUTPUT_DIR}/input/admin-boundaries/zip/gadm36_gpkg.zip"
    shell:
        """
        wget -O {output} https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_gpkg.zip
        """

rule unzip_gadm:
    input:
        "{OUTPUT_DIR}/input/admin-boundaries/zip/gadm36_gpkg.zip"
    output:
        "{OUTPUT_DIR}/input/admin-boundaries/gadm36.gpkg"
    params:
        dirname=lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        unzip -o {input} -d {params.dirname}
        """

# Not very DRY, but okay for now
rule download_ne_50m:
    output:
        "{OUTPUT_DIR}/input/admin-boundaries/zip/ne_50m.zip"
    shell:
        """
        wget -O {output} https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/cultural/ne_50m_admin_0_countries.zip
        """

rule unzip_ne_50m:
    input:
        "{OUTPUT_DIR}/input/admin-boundaries/zip/ne_50m.zip"
    output:
        directory("{OUTPUT_DIR}/input/admin-boundaries/ne_50m/")
    shell:
        """
        unzip -o {input} -d {output}
        """

"""
Test with:
snakemake --cores 1 results/input/admin-boundaries/gadm36.gpkg
"""

# download admin boundaries per country
rule download_gadm_by_country:
    output:
        os.path.join(config['output_dir'], "input", "adminboundaries", "gadm36_{code}.gpkg"),
    shell:
        """
        wget https://biogeo.ucdavis.edu/data/gadm3.6/gpkg/gadm36_{wildcards.code}_gpkg.zip \
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
        wget https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_levels_gpkg.zip \
            --directory-prefix={config['output_dir']}/input/adminboundaries
        unzip -o {config['output_dir']}/input/adminboundaries/gadm36_levels_gpkg.zip -d {config['output_dir']}/input/adminboundaries
        """
