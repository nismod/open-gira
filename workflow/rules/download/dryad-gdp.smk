"""Download Dryad gridded GDP

Reference
---------
https://doi.org/10.5061/dryad.dk1j0
"""

GDP_DATASETS = expand(
    os.path.join(config["output_dir"], "input", "GDP", "{filename}"),
    filename=[
        "admin_areas_GDP_HDI.nc",
        "GDP_per_capita_PPP_1990_2015_v2.nc",
        "GDP_PPP_30arcsec_v3.nc",
        "GDP_PPP_1990_2015_5arcmin_v2.nc",
        "HDI_1990_2015_v2.nc",
        "pedigree_GDP_per_capita_PPP_1990_2015_v2.nc",
        "pedigree_HDI_1990_2015_v2.nc",
    ],
)


rule download_GDP:
    output:
        GDP_DATASETS,
    shell:
        f"""
        mkdir -p {config['output_dir']}/input/GDP
        cd {config['output_dir']}/input/GDP
        wget https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.dk1j0/download \
            --content-disposition
        unzip -o doi_10.5061_dryad.dk1j0__v2.zip
        """
