"""
Convert .geoparquet files into raster files by averaging road lengths by
cell index. Cell indices are then used to build the raster file.
"""


rule make_exposure_img:
    conda: "../../../../environment.yml"
    input:
        hazard_dir="{OUTPUT_DIR}/exposure/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/raster/",
        coastline="{OUTPUT_DIR}/input/coastlines/ne_10m_ocean/",
        boundaries="{OUTPUT_DIR}/input/admin-boundaries/ne_50m/",
    output:
        directory("{OUTPUT_DIR}/exposure/{DATASET}_{FILTER_SLUG}/{HAZARD_SLUG}/img"),
    params:
        # Put the extensions in params because otherwise SNAKEMAKE doesn't find them
        coastline="ne_10m_ocean.shp",
        boundaries="ne_50m_admin_0_countries.shp",
    script:
        "../../../scripts/make_exposure_img.py"


"""
Test with:
snakemake --cores all results/exposure/tanzania-mini_filter-road/hazard-aqueduct-river/img/
"""
