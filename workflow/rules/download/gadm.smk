"""Download GADM boundaries

Reference
---------
https://gadm.org/data.html
"""

out_adminboundaries = os.path.join(DATA_DIR, "adminboundaries", "gadm36.gpkg")
out_adminboundaries_codes = expand(
    os.path.join(DATA_DIR, "adminboundaries", "gadm36_{code}.gpkg"), code=COUNTRY_CODES
)


rule download_gadm:
    output:
        out_adminboundaries,
    shell:
        """
        wget https://biogeo.ucdavis.edu/data/gadm3.6/gadm36_gpkg.zip \
            --directory-prefix=data/adminboundaries
        unzip -o gadm36_gpkg.zip
        """


# download admin boundaries per country
rule download_gadm_by_country:
    output:
        os.path.join(DATA_DIR, "adminboundaries", "gadm36_{code}.gpkg"),
    shell:
        """
        wget https://biogeo.ucdavis.edu/data/gadm3.6/gpkg/gadm36_{wildcards.code}_gpkg.zip \
            --directory-prefix=data/adminboundaries
        unzip -o data/adminboundaries/gadm36_{wildcards.code}_gpkg.zip -d data/adminboundaries
        """


# rule download_all_per_country__masterrule:
#     input:
#         expand(os.path.join(DATA_DIR, "adminboundaries", "gadm36_{code}.gpkg"), code=COUNTRY_CODES)