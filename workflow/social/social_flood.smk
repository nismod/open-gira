rule clip_flood_data:
    """
    This rule is adjusted from the "trim_rasters" rule in "trim_hazard_data.smk"
    This is a temporary rule, can get it working with trim rasters eventually. 
    TODO: delete rule and get workflow working with trim_rasters.
    """
    input:
        flood_file="{OUTPUT_DIR}/input/{HAZARD_SLUG}/merged/jrc_global_flood_RP{RP}.tif",
        json_file="{OUTPUT_DIR}/json/{DATASET}.json",
    wildcard_constraints:
        RP="10|20|50|75|100|200|500"
    output:
        trimmed_flood_file="{OUTPUT_DIR}/input/{HAZARD_SLUG}/trimmed/{DATASET}/jrc_global_flood_RP{RP}.tif",
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_flood_file})

        # pull out bounding box coords into bash array
        declare -a "COORDS=($(jq '.extracts[0].bbox | @sh' {input.json_file}))"

        # now trim the raster

        # a plain gdalwarp call will trim the dataset successfully
        # however, the output will be tens of times larger per unit area than the original
        # using gdal_translate we can trim _and_ compress
        # https://trac.osgeo.org/gdal/wiki/UserDocs/GdalWarp#GeoTIFFoutput-coCOMPRESSisbroken

        # first, generate a VRT format job specification (a short XML file) with gdalwarp
        JOB_SPEC=$(mktemp)  # file in /tmp
        gdalwarp -te ${{COORDS[@]}} -of VRT {input.flood_file} $JOB_SPEC

        # then use gdal_translate to execute the job as specified
        gdal_translate -co compress=lzw $JOB_SPEC {output.trimmed_flood_file}

        # clean up job specification file
        rm $JOB_SPEC
        """
"""
Test with
snakemake -c1 results/input/hazard-jrc-river/trimmed/kenya-latest/jrc_global_flood_RP10.tif
"""

rule flood_exposed_rwi:
    """
    This rule overlays the flood depth rasters with the rwi rasters and returnes the rwi
    value wherever flood depth > 1.
    """
    input:
        flood_file="{OUTPUT_DIR}/input/{HAZARD_SLUG}/trimmed/{DATASET}/jrc_global_flood_RP{RP}.tif",
        rwi_file="{OUTPUT_DIR}/input/rwi/{DATASET}/rwi.tif",
    wildcard_constraints:
        RP="10|20|50|75|100|200|500"
    output:
        rwi_exposure="{OUTPUT_DIR}/social/rwi/gridded_flood_exposure/{HAZARD_SLUG}/{DATASET}/rwi_exposure_RP{RP}.tif"
    script:
        "./flood_exposed_social.py"
"""
Test with
snakemake -c1 results/social/rwi/gridded_flood_exposure/hazard-jrc-river/kenya-latest/rwi_exposure_RP10.tif
"""   

