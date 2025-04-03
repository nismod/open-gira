"""Download Dryad gridded GDP

See API docs at https://datadryad.org/api

- dataset record https://datadryad.org/api/v2/datasets/doi:10.5061%2Fdryad.dk1j0
  - find latest version id
- version files https://datadryad.org/api/v2/versions/52523/files
  - find list of files with links and ids

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
        output_dir=$(dirname {output.admin})

        for file_id in 241946 241947 241948 241949 241950 241951 241953 241958
        do
            wget https://datadryad.org/api/v2/files/$file_id/download \
                -nc \
                --content-disposition \
                --directory-prefix=$output_dir
        """

"""
Test with:
snakemake -c1 -- results/input/GDP/GDP_per_capita_PPP_1990_2015_v2.nc
"""
