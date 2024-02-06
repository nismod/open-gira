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
