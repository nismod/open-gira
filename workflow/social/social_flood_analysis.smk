"""
Standalone snakemake containing all workflow for Mark's paper idea
"""

rule clip_rwi:
    """
    Clip rwi raster to analysis region
    """
    input:
        raw_rwi_file="{OUTPUT_DIR}/input/rwi/rwi.tif",
        json_file="{OUTPUT_DIR}/json/{DATASET}.json",
    output:
        trimmed_rwi_file="{OUTPUT_DIR}/mark_paper/rwi/{DATASET}/rwi.tif",
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_rwi_file})
        
        # pull out bounding box coords into bash array
        declare -a "COORDS=($(jq '.extracts[0].bbox | @sh' {input.json_file}))"

        # now trim the raster

        # Reproject rwi dataset (to EPSG:4326)
        gdalwarp \
            -te ${{COORDS[@]}} \
            -t_srs EPSG:4326 \
            -r nearest \
            -dstnodata -999 \
            {input.raw_rwi_file} \
            {output.trimmed_rwi_file}
        """
"""
Test with
snakemake -c1 results/mark_paper/rwi/kenya-latest/rwi.tif
"""

rule rwi_flood_risk:
    """
    This rule overlays the flood hazard map with the rwi map and extracts
    flood exposed rwi grid cells across various vulnerability thresholds (according to damage curve).
    Analysis is conducted at resolution of the rwi data (2.4 km). Flood data is 
    resampled accordingly.
    """
    input:
        flood_file="{OUTPUT_DIR}/input/{HAZARD_SLUG}/trimmed/{DATASET}/jrc_global_flood_RP{RP}.tif",
        rwi_file="{OUTPUT_DIR}/mark_paper/rwi/{DATASET}/rwi.tif",
    wildcard_constraints:
        RP="10|20|50|75|100|200|500"
    output:
        rwi_vhigh_risk="{OUTPUT_DIR}/mark_paper/rwi_risk/{HAZARD_SLUG}/{DATASET}/very_high/rwi_risk_RP{RP}.tif",
        rwi_high_risk="{OUTPUT_DIR}/mark_paper/rwi_risk/{HAZARD_SLUG}/{DATASET}/high/rwi_risk_RP{RP}.tif",
        rwi_medium_risk="{OUTPUT_DIR}/mark_paper/rwi_risk/{HAZARD_SLUG}/{DATASET}/medium/rwi_risk_RP{RP}.tif",
        rwi_low_risk="{OUTPUT_DIR}/mark_paper/rwi_risk/{HAZARD_SLUG}/{DATASET}/low/rwi_risk_RP{RP}.tif",
        rwi_exposure="{OUTPUT_DIR}/mark_paper/rwi_risk/{HAZARD_SLUG}/{DATASET}/exposure/rwi_exposure_RP{RP}.tif",
    script:
        "./social_flood_risk.py"
"""
Test with
snakemake -c1 results/mark_paper/rwi_risk/hazard-jrc-river/kenya-latest/very_high/rwi_risk_RP10.tif
"""

rule aggregate_population_data:
    """
    This rule implements area-weighted aggregation to resample population data to the resolution 
    of the RWI, while preserving population estimate accuracy as much as possible.
    """
    input:
        pop_file="{OUTPUT_DIR}/input/ghsl/{DATASET}/GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif",
        rwi_file="{OUTPUT_DIR}/mark_paper/rwi/{DATASET}/rwi.tif",
    output:
        agg_pop_file="{OUTPUT_DIR}/mark_paper/pop/aggregated/{DATASET}/pop.tif"
    script:
        "./aggregate_population.py"
"""
Test with
snakemake -c1 results/mark_paper/pop/aggregated/kenya-latest/pop.tif"
"""