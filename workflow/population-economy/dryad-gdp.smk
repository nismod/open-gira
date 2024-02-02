"""Download Dryad gridded GDP

Reference
---------
https://doi.org/10.5061/dryad.dk1j0
"""

rule download_GDP:
    output:
        admin = "{OUTPUT_DIR}/input/GDP/admin_areas_GDP_HDI.nc",
        gdp_pc = "{OUTPUT_DIR}/input/GDP/GDP_per_capita_PPP_1990_2015_v2.nc",
        gdp_ppp_30arcsec = "{OUTPUT_DIR}/input/GDP/GDP_PPP_30arcsec_v3.nc",
        gdp_ppp_5arcmin = "{OUTPUT_DIR}/input/GDP/GDP_PPP_1990_2015_5arcmin_v2.nc",
        gdp_pedigree_pc = "{OUTPUT_DIR}/input/GDP/pedigree_GDP_per_capita_PPP_1990_2015_v2.nc",
        hdi = "{OUTPUT_DIR}/input/GDP/HDI_1990_2015_v2.nc",
        hdi_pedigree = "{OUTPUT_DIR}/input/GDP/pedigree_HDI_1990_2015_v2.nc",
    shell:
        """
        mkdir -p {wildcards.OUTPUT_DIR}/input/GDP
        cd {wildcards.OUTPUT_DIR}/input/GDP
        wget https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.dk1j0/download \
            --content-disposition
        unzip -o doi_10.5061_dryad.dk1j0__v2.zip
        """

"""
Test with:
snakemake -c1 -- results/input/GDP/GDP_per_capita_PPP_1990_2015_v2.nc
"""
