"""
Convert .geoparquet files into raster files by averaging road lengths by
cell index. Cell indices are then used to build the raster file.
"""

# https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#data-dependent-conditional-execution
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.make_exposure_img.get(**wildcards).output[0]
    tifs = glob_wildcards(os.path.join(checkpoint_output, "{tif}.tif"))
    input = expand(os.path.join(checkpoint_output, "{tif}.tif"), tif=tifs.tif)
    return input

rule make_exposure_img:
    input:
        geoparquet="{OUTPUT_DIR}/{DATASET}_{FILTER_SLUG}_{HAZARD_SLUG}.geoparquet",
        hazard=aggregate_input,
        coastline="{OUTPUT_DIR}/input/ne_10m_ocean/ne_10m_ocean.shp",
        boundaries="{OUTPUT_DIR}/input/admin-boundaries/ne_50m/ne_50m_admin_0_countries.shp"
    output:
        directory("{OUTPUT_DIR}/exposure/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/img")
    script:
        "../scripts/make_exposure_img.py"

"""
Test with:
snakemake --cores all results/exposure/tanzania-mini_filter-highway-core/hazard-aqueduct-river/img/
"""
