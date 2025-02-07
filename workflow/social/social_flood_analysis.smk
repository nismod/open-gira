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