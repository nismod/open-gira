"""Download GADM boundaries

Reference
---------
https://gadm.org/data.html
"""

# TODO follow config data_dir in scripts
rule download_gadm:
    output:
        os.path.join(config['data_dir'], "adminboundaries", "gadm36.gpkg"),
    shell:
        """
        wget https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_gpkg.zip \
            --directory-prefix=data/adminboundaries
        unzip -o gadm36_gpkg.zip
        """


rule download_gadm_levels:
    output:
        os.path.join(config['data_dir'], "adminboundaries", "gadm36_levels.gpkg"),
    shell:
        """
        wget https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_levels_gpkg.zip \
            --directory-prefix=data/adminboundaries
        unzip -o gadm36_levels_gpkg.zip
        """


rule download_gadm_by_country:
    output:
        os.path.join(config['data_dir'], "adminboundaries", "gadm36_{code}.gpkg"),
    shell:
        """
        wget https://biogeo.ucdavis.edu/data/gadm3.6/gpkg/gadm36_{wildcards.code}_gpkg.zip \
            --directory-prefix=data/adminboundaries
        unzip -o data/adminboundaries/gadm36_{wildcards.code}_gpkg.zip \
            -d data/adminboundaries
        """