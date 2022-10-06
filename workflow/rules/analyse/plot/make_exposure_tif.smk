"""
Convert .geoparquet files into raster files by averaging road lengths by
cell index. Cell indices are then used to build the raster file.
"""

# https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#data-dependent-conditional-execution
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.trim_hazard_data.get(**wildcards).output[0]
    tifs = glob_wildcards(os.path.join(checkpoint_output, "{tif}.tif"))
    input = expand(os.path.join(checkpoint_output, "{tif}.tif"), tif=tifs.tif)
    return input


rule make_exposure_tif:
    input:
        geoparquet="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}_{HAZARD_SLUG}.geoparquet",
        hazard=aggregate_input,
    output:
        directory("{OUTPUT_DIR}/exposure/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/raster/"),
    script:
        "../../../scripts/make_exposure_tif.py"


"""
Test with:
snakemake --cores all results/exposure/tanzania-mini_filter-road/hazard-aqueduct-river/raster/
"""
