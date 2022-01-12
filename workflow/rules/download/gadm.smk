"""Download GADM boundaries

Reference
---------
https://gadm.org/data.html
"""

out_adminboundaries = os.path.join(DATA_DIR, "adminboundaries", "gadm36.gpkg")
out_adminboundaries_levels = os.path.join(DATA_DIR, "adminboundaries", "gadm36_levels.gpkg")

rule download_gadm:
    output:
        out_adminboundaries,
    shell:
        """
        wget https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_gpkg.zip \
            --directory-prefix=data/adminboundaries
        unzip -o gadm36_gpkg.zip
        """


rule download_gadm_levels:
    output:
        out_adminboundaries_levels,
    shell:
        """
        wget https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_levels_gpkg.zip \
            --directory-prefix=data/adminboundaries
        unzip -o gadm36_levels_gpkg.zip
        """
