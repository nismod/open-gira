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

rule resample_flood_data:
    """
    This rule resamples flood data to the resolution of the RWI datasets.
    Here resample by preserving the maximum flooded cell within each gridded cell.
    """
    input:
        flood_file="{OUTPUT_DIR}/input/{HAZARD_SLUG}/trimmed/{DATASET}/jrc_global_flood_RP{RP}.tif",
        rwi_file="{OUTPUT_DIR}/mark_paper/rwi/{DATASET}/rwi.tif"
    wildcard_constraints:
        RP="10|20|50|75|100|200|500"
    output:
        resampled_flood_file = "{OUTPUT_DIR}/input/{HAZARD_SLUG}/rwi_resampled/{DATASET}/jrc_global_flood_RP{RP}.tif"
    script:
        "./resample_flood_data.py"
"""
Test with
snakemake -c1 results/input/hazard-jrc-river/rwi_resampled/kenya-latest/jrc_global_flood_RP10.tif
"""

rule relative_risk_flood:
    """
    This rule calculates the relative risk per grid cell according to a vulnerability curve.
    The two curves that we use are either the JRC residential (average of Asia, Africa, and SAmerica)
    or the Bernhofen curve (refugee paper)
    """
    input:
        flood_file="{OUTPUT_DIR}/input/{HAZARD_SLUG}/rwi_resampled/{DATASET}/jrc_global_flood_RP{RP}.tif"
    output:
        risk_file="{OUTPUT_DIR}/mark_paper/relative_risk/{HAZARD_SLUG}/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP{RP}.tif"
    wildcard_constraints:
        VULN_CURVE="BERN|JRC",
        RP="10|20|50|75|100|200|500"
    script:
        "./relative_flood_risk.py"
""" 
Test with
snakemake -c1 results/mark_paper/relative_risk/hazard-jrc-river/kenya-latest/JRC/jrc_global_flood_RP10.tif
"""

rule average_annual_risk_jrc:
    """
    This rule calculates one layer (average annual relative risk) given a set of 
    return period flood maps. This rule is written specifically for the JRC maps.
    Also calcualtes the "protected" AAR - by integrating flopros layer
    """
    input:
        flood_rp_10="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP10.tif",
        flood_rp_20="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP20.tif",
        flood_rp_50="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP50.tif",
        flood_rp_75="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP75.tif",
        flood_rp_100="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP100.tif",
        flood_rp_200="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP200.tif",
        flood_rp_500="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP500.tif",
    output:
        flood_aar="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_AAR.tif",
    wildcard_constraints:
        VULN_CURVE="BERN|JRC",
    script:
        "./average_annual_risk.py"
"""
Test with
snakemake -c1 results/mark_paper/relative_risk/hazard-jrc-river/kenya-latest/JRC/jrc_global_flood_AAR.tif
"""

rule download_ghsmod:
    """
    Download the GHS-MOD dataset from GHSL website.
    Will be used to classify urban and rural areas.
    Also resampling to WGS84 (30 arcsecs) using nearest neighbor
    """
    output:
        "{OUTPUT_DIR}/input/ghsl/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V2_0.tif"
    shell:
        """
        output_dir=$(dirname {output})

        mkdir -p $output_dir

        wget -nc https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_SMOD_GLOBE_R2023A/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000/V2-0/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V2_0.zip \
            --directory-prefix=$output_dir

        unzip -o $output_dir/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V2_0.zip \
            -d $output_dir
        """

rule clip_and_resample_GHSMOD:
    """
    Clips the GHSMOD dataset, reprojects (to WGS84), and resamples to RWI resolution 
    using the nearest neighbor approach
    """
    input:
        ghs_mod="{OUTPUT_DIR}/input/ghsl/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V2_0.tif",
        json_file="{OUTPUT_DIR}/json/{DATASET}.json",
        rwi_file="{OUTPUT_DIR}/mark_paper/rwi/{DATASET}/rwi.tif"
    output:
        ghs_mod_clipped = "{OUTPUT_DIR}/mark_paper/ghsl/{DATASET}/GHS_SMOD_E2020.tif"
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.ghs_mod_clipped})
        
        # pull out bounding box coords into bash array
        declare -a "COORDS=($(jq '.extracts[0].bbox | @sh' {input.json_file}))"

        # Pull reference pixel resolution
        tr=$(gdalinfo {input.rwi_file} | grep "Pixel Size" | head -n1 | sed -e 's/.*(//' -e 's/).*//' | awk -F, '{{print $1, $2}}')

        # now trim the raster

        # Reproject rwi dataset (to EPSG:4326)
        gdalwarp \
            -te ${{COORDS[@]}} \
            -t_srs EPSG:4326 \
            -r nearest \
            -tr $tr \
            {input.ghs_mod} \
            {output.ghs_mod_clipped}
        """
"""
Test with 
snakemake -c1 results/mark_paper/ghsl/kenya-latest/GHS_SMOD_E2020.tif
""" 

rule download_flopros:
    """
    Download and extract the FLOPROS dataset.
    This rule outputs the shapefile needed for further processing.
    """
    output:
        shp="{OUTPUT_DIR}/input/flopros/Scussolini_etal_Suppl_info/FLOPROS_shp_V1/FLOPROS_shp_V1.shp"
    shell:
        """
        output_dir=$(dirname $(dirname $(dirname {output.shp})))
        mkdir -p $output_dir
        wget -nc http://dx.doi.org/10.5194/nhess-16-1049-2016-supplement \
            --directory-prefix=$output_dir
        unzip -o $output_dir/nhess-16-1049-2016-supplement
            -d $output_dir
        """
"""
Test with 
snakemake -c1 results/input/flopros/Scussolini_etal_Suppl_info/FLOPROS_shp_V1/FLOPROS_shp_V1.shp"
""" 

rule trim_flopros:
    """
    Trim FLOPROS dataset to region of interest
    """
    input:
        flopros = "{OUTPUT_DIR}/input/flopros/Scussolini_etal_Suppl_info/FLOPROS_shp_V1/FLOPROS_shp_V1.shp",
        json_file="{OUTPUT_DIR}/json/{DATASET}.json",
    output:
        flopros_clipped = "{OUTPUT_DIR}/mark_paper/flopros/{DATASET}/flopros.shp"
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.flopros_clipped})

        # pull out bounding box coords into bash array
        declare -a "COORDS=($(jq '.extracts[0].bbox | @sh' {input.json_file}))"
        
        ogr2ogr -f "ESRI Shapefile" {output.flopros_clipped} {input.flopros} -clipsrc ${{COORDS[@]}}
        """
"""
Test with 
snakemake -c1 results/mark_paper/flopros/kenya-latest/flopros.shp
"""

rule rasterize_flopros:
    """
    We are going to rasterize the flopros layer (using the merged protection value column).
    """
    input:
        flopros="{OUTPUT_DIR}/mark_paper/flopros/{DATASET}/flopros.shp",
        rwi_file="{OUTPUT_DIR}/mark_paper/rwi/{DATASET}/rwi.tif"
    output:
        rasterized_flopros="{OUTPUT_DIR}/mark_paper/flopros/{DATASET}/flopros.tif"
    script:
        "./rasterize_flopros.py"
"""
Test with 
snakemake -c1 results/mark_paper/flopros/kenya-latest/flopros.tif
""" 

rule average_annual_risk_jrc_protected:
    """
    This rule calculates one layer (average annual relative risk) given a set of 
    return period flood maps. This rule is written specifically for the JRC maps.
    Also calcualtes the "protected" AAR - by integrating flopros layer
    """
    input:
        flood_rp_10="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP10.tif",
        flood_rp_20="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP20.tif",
        flood_rp_50="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP50.tif",
        flood_rp_75="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP75.tif",
        flood_rp_100="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP100.tif",
        flood_rp_200="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP200.tif",
        flood_rp_500="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP500.tif",
        flopros="{OUTPUT_DIR}/mark_paper/flopros/{DATASET}/flopros.tif",
    output:
        flood_aar="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_AAR_protected.tif",
    wildcard_constraints:
        VULN_CURVE="BERN|JRC",
    script:
        "./average_annual_risk_protected.py"
"""
Test with
snakemake -c1 results/mark_paper/relative_risk/hazard-jrc-river/kenya-latest/JRC/jrc_global_flood_AAR_protected.tif
"""

rule concentration_index:
    """
    This rule calcualtes flood risk concentration indices for different admin regions.
    Here calculating it using the JRC damage curve...
    """
    input:
        admin_areas = "{OUTPUT_DIR}/input/admin-boundaries/{ADMIN_SLUG}.geoparquet",
        json_file="{OUTPUT_DIR}/json/{DATASET}.json",
        rwi_file="{OUTPUT_DIR}/mark_paper/rwi/{DATASET}/rwi.tif",
        pop_file="{OUTPUT_DIR}/mark_paper/pop/aggregated/{DATASET}/pop.tif",
        urban_file="{OUTPUT_DIR}/mark_paper/ghsl/{DATASET}/GHS_SMOD_E2020.tif",
        risk_file="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/JRC/jrc_global_flood_AAR.tif",
        risk_protected_file="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/JRC/jrc_global_flood_AAR_protected.tif"
    output:
        regional_CI = "{OUTPUT_DIR}/mark_paper/concentration_index/{DATASET}/concentration_index_{ADMIN_SLUG}.gpkg",
    script:
        "./concentration_index.py"
"""
Test with
snakemake -c1 results/mark_paper/concentration_index/kenya-latest/concentration_index_admin-level-1.gpkg
"""

rule plot_concentration_curve:
    """
    This rule plots the concentration curve for the chosen admin region.
    Here calculating it using the JRC damage curve...
    """
    input:
        admin_areas = "{OUTPUT_DIR}/input/admin-boundaries/{ADMIN_SLUG}.geoparquet",
        json_file="{OUTPUT_DIR}/json/{DATASET}.json",
        rwi_file="{OUTPUT_DIR}/mark_paper/rwi/{DATASET}/rwi.tif",
        pop_file="{OUTPUT_DIR}/mark_paper/pop/aggregated/{DATASET}/pop.tif",
        urban_file="{OUTPUT_DIR}/mark_paper/ghsl/{DATASET}/GHS_SMOD_E2020.tif",
        risk_file="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/JRC/jrc_global_flood_AAR.tif",
        risk_protected_file="{OUTPUT_DIR}/mark_paper/relative_risk/hazard-jrc-river/{DATASET}/JRC/jrc_global_flood_AAR_protected.tif"
    output:
        figure_directory = directory("{OUTPUT_DIR}/mark_paper/concentration_index/{DATASET}/plots/{ADMIN_SLUG}/"),
    script:
        "./plot_concentration_curves.py"
"""
Test with
snakemake -c1 results/mark_paper/concentration_index/kenya-latest/plots/admin-level-0
"""