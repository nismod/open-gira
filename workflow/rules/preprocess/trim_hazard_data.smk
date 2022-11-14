# chop hazard data according to overall bounding box (stored in dataset's JSON file)

from json import load
from pathlib import Path


# https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#data-dependent-conditional-execution
def generate_trimmed_tif_paths(wildcards):
    # where the raw (untrimmed) tiffs live
    raw_folder = Path(checkpoints.download_hazard_datasets.get(**wildcards).output.raw_folder)
    # where the trimmed tiffs for a given DATASET go
    dataset_folder: Path = raw_folder.parent / wildcards.DATASET
    # file basenames to create (trimmed) full paths for
    tiff_basenames, = glob_wildcards(raw_folder / "{BASENAME}.tif")
    # full paths for trimmed tiff files
    return expand(dataset_folder / "{basename}.tif", basename=tiff_basenames)


rule trim_hazard_data:
    """
    Collect all of the trimmed raster filenames pertinent to a given dataset

    Essentially a target rule
    """
    input:
        trimmed_rasters=generate_trimmed_tif_paths,


rule trim_raster:
    """
    Trim a single raster according to a JSON specification
    """
    input:
        raw_raster_file="{OUTPUT_DIR}/input/{HAZARD_SLUG}/raw/{TIFF_FILE}",
        json_file="{OUTPUT_DIR}/json/{DATASET}.json",
    output:
        trimmed_raster_file="{OUTPUT_DIR}/input/{HAZARD_SLUG}/{DATASET}/{TIFF_FILE}",
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_raster_file})

        # pull out bounding box coords into bash array
        declare -a "COORDS=($(jq '.extracts[0].bbox | @sh' {input.json_file}))"

        # now trim the raster

        # a plain gdalwarp call will trim the dataset successfully
        # however, the output will be tens of times larger per unit area than the original
        # using gdal_translate we can trim _and_ compress
        # https://trac.osgeo.org/gdal/wiki/UserDocs/GdalWarp#GeoTIFFoutput-coCOMPRESSisbroken

        # first, generate a VRT format job specification (a short XML file) with gdalwarp
        JOB_SPEC=$(mktemp)  # file in /tmp
        gdalwarp -te ${{COORDS[@]}} -of VRT {input.raw_raster_file} $JOB_SPEC

        # then use gdal_translate to execute the job as specified
        gdal_translate -co compress=lzw $JOB_SPEC {output.trimmed_raster_file}

        # clean up job specification file
        rm $JOB_SPEC
        """

"""
Test with:
snakemake --cores 1 results/input/hazard-aqueduct-river/planet-latest/inunriver_rcp8p5_00000NorESM1-M_2080_rp00005.tif
"""
