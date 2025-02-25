"""
Downloading land cover/usage maps
"""

rule download_glob_cover_2009:
    """
    Download file from ESA remote.
    """
    output:
        zip_file = "{OUTPUT_DIR}/input/land_cover/Globcover2009_V2.3_Global_.zip"
    shell:
        """
        wget http://due.esrin.esa.int/files/Globcover2009_V2.3_Global_.zip -O {output.zip_file}
        """

rule unzip_glob_cover_2009:
    """
    Extract contents of archive.

    Data file is named: GLOBCOVER_L4_200901_200912_V2.3.tif, in output.unzip_dir.
    """
    input:
        zip_file = rules.download_glob_cover_2009.output.zip_file
    output:
        unzip_dir = directory("{OUTPUT_DIR}/input/land_cover/glob_cover_2009"),
        raster = "{OUTPUT_DIR}/input/land_cover/glob_cover_2009/GLOBCOVER_L4_200901_200912_V2.3.tif"
    shell:
        """
        unzip {input.zip_file} -d {output.unzip_dir}
        """

rule download_copernicus_global_dynamic_land_cover:
    """Copernicus Global Land Service: Land Cover 100m

    Near real time epoch 2019 from the Collection 3 of annual, global 100m land cover maps.

    Other available epochs: 2015  2016  2017  2018

    Produced by the global component of the Copernicus Land Service, derived from PROBA-V satellite observations and ancillary datasets.

    The maps include:
    - a main discrete classification with 23 classes aligned with UN-FAO's Land Cover Classification System,
    - a set of versatile cover fractions: percentage (%) of ground cover for the 10 main classes
    - a forest type layer
    - quality layers on input data density and on the confidence of the detected land cover change


    Marcel Buchhorn, Bruno Smets, Luc Bertels, Bert De Roo, Myroslava Lesiv,
    Nandin-Erdene Tsendbazar, Martin Herold, & Steffen Fritz. (2020). Copernicus
    Global Land Service: Land Cover 100m: collection 3: epoch 2019: Globe
    (V3.0.1) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.3939050
    """
    output:
        links="{OUTPUT_DIR}/input/land_cover/copernicus_global_dynamic_land_cover/links.txt",
        crop_tif="{OUTPUT_DIR}/input/land_cover/copernicus_global_dynamic_land_cover/PROBAV_LC100_global_v3.0.1_2019-nrt_Crops-CoverFraction-layer_EPSG-4326.tif",
        tree_tif="{OUTPUT_DIR}/input/land_cover/copernicus_global_dynamic_land_cover/PROBAV_LC100_global_v3.0.1_2019-nrt_Tree-CoverFraction-layer_EPSG-4326.tif",
        built_tif="{OUTPUT_DIR}/input/land_cover/copernicus_global_dynamic_land_cover/PROBAV_LC100_global_v3.0.1_2019-nrt_BuiltUp-CoverFraction-layer_EPSG-4326.tif",
    shell:
        """
        pushd $(dirname {output.links})
            zenodo_get -w links.txt --record=3939050
            wget -nc -i links.txt
            md5sum -c md5sums.txt
        popd
        """
