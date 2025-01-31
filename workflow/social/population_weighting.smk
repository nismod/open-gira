rule download_pop:
    """
    Note: this is just the "download_ghsl" rule from "ghsl-pop.smk"
    I have re-used it here to adjust to using the WGS84 population maps
    to avoid messing with too many rules on this branch. Can adjust (remove)
    this rule down the line.
    """
    output:
        "{OUTPUT_DIR}/input/ghsl/GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif"
    shell:
        """
        output_dir=$(dirname {output})

        mkdir -p $output_dir

        wget -nc https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E2020_GLOBE_R2023A_4326_30ss/V1-0/GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.zip \
            --directory-prefix=$output_dir

        unzip -o $output_dir/GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.zip \
            -d $output_dir
        """

rule clip_rasters:
    """
    Clip population and rwi rasters for population weighting. This rule is adjusted from the 
    "trim_rasters" rule in "trim_hazard_data.smk"
    """
    input:
        raw_rwi_file="{OUTPUT_DIR}/input/rwi/rwi.tif",
        raw_pop_file="{OUTPUT_DIR}/input/ghsl/GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif",
        json_file="{OUTPUT_DIR}/json/{DATASET}.json",
    output:
        trimmed_rwi_file="{OUTPUT_DIR}/input/rwi/{DATASET}/rwi.tif",
        trimmed_pop_file="{OUTPUT_DIR}/input/ghsl/{DATASET}/GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif",
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_rwi_file})
        
        # pull out bounding box coords into bash array
        declare -a "COORDS=($(jq '.extracts[0].bbox | @sh' {input.json_file}))"

        # now trim the raster

        # calling gdalwarp (note that trim_rasters original rule uses gdal_translate due to file sizes)
        # First for pop dataset
        gdalwarp -te ${{COORDS[@]}} {input.raw_pop_file} {output.trimmed_pop_file}
        # For rwi dataset we will also be reprojecting (to EPSG:4326) and resampling (to 0.00833 deg ~ 1km) to match pop data
        gdalwarp \
            -te ${{COORDS[@]}} \
            -t_srs EPSG:4326 \
            -r bilinear \
            -dstnodata -999 \
            -tr 0.008333 0.008333 \
            {input.raw_rwi_file} \
            {output.trimmed_rwi_file}
        """

rule population_weight_rwi_at_admin_area:
    input:
        admin_areas = "{OUTPUT_DIR}/input/admin-boundaries/{ADMIN_SLUG}.geoparquet",
        rwi_file = "{OUTPUT_DIR}/input/rwi/{DATASET}/rwi.tif",
        pop_file = "{OUTPUT_DIR}/input/ghsl/{DATASET}/GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif",
        bbox_file = "{OUTPUT_DIR}/json/{DATASET}.json",
    output:
        pop_weighted_rwi = "{OUTPUT_DIR}/social/rwi/population_weighted_adm/{DATASET}/pop_weighted_rwi_{ADMIN_SLUG}.geoparquet",
    script:
        "./population_weight_at_admin_area.py"

"""
Test with
snakemake -c1 results/social/rwi/population_weighted_adm/kenya-latest/pop_weighted_rwi_admin-level-1.geoparquet
"""