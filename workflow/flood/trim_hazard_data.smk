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
        tiff="{OUTPUT_DIR}/input/{HAZARD_SLUG}/raw/{TIFF_FILE}",
        json="{OUTPUT_DIR}/json/{DATASET}.json",
    output:
        tiff="{OUTPUT_DIR}/input/{HAZARD_SLUG}/{DATASET}/{TIFF_FILE}",
    shell:
        """
        set -ex

        mkdir -p $(dirname {output.tiff})

        # pull out bounding box coords into bash array
        # bbox is [xmin, ymin, xmax, ymax]
        declare -a "COORDS=($(jq '.extracts[0].bbox | @sh' {input.json}  | sed 's/"//g'))"

        # then use gdal_translate with projwin to avoid warping
        # projwin is [ulx, uly, lrx, lry]
        gdal_translate \
            -projwin ${{COORDS[0]}} ${{COORDS[3]}} ${{COORDS[2]}} ${{COORDS[1]}} \
            -co COMPRESS=LZW \
            -co BIGTIFF=IF_SAFER \
            --config GDAL_CACHEMAX 50% \
            {input.tiff} {output.tiff}

        """

"""
Test with:
snakemake --cores 1 results/input/hazard-aqueduct-river/planet-latest/inunriver_rcp8p5_00000NorESM1-M_2080_rp00005.tif
"""
