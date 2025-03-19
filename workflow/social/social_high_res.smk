"""
Standalone snakemake for conducting Mark's paper analysis at high (3as) resolution
"""

rule clip_3as_pop_raster:
    """
    Clip 3 as population raster. This rule is adjusted from the 
    "trim_rasters" rule in "trim_hazard_data.smk"
    """
    input:
        raw_pop_file="{OUTPUT_DIR}/input/ghsl/GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif",
        json_file="{OUTPUT_DIR}/json/{DATASET}.json",
    output:
        trimmed_pop_file="{OUTPUT_DIR}/mark_paper/high_resolution/ghsl/{DATASET}/GHS_POP_E2020.tif",
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_pop_file})
        
        # pull out bounding box coords into bash array
        declare -a "COORDS=($(jq '.extracts[0].bbox | @sh' {input.json_file}))"

        # now trim the raster

        # calling gdalwarp (note that trim_rasters original rule uses gdal_translate due to file sizes)
        # First for pop dataset
        gdalwarp -te ${{COORDS[@]}} {input.raw_pop_file} {output.trimmed_pop_file}
        """
""" 
Test with
snakemake -c1 results/mark_paper/high_resolution/ghsl/kenya-latest/GHS_POP_E2020.tif
"""

rule resample_rwi_and_ghsmod:
    """
    Resample the rwi and ghsmod data to the resolution of the population data
    """
    input:
        pop_file="{OUTPUT_DIR}/mark_paper/high_resolution/ghsl/{DATASET}/GHS_POP_E2020.tif",
        rwi_file="{OUTPUT_DIR}/mark_paper/rwi/{DATASET}/rwi.tif",
        json_file="{OUTPUT_DIR}/json/{DATASET}.json",
        ghs_mod="{OUTPUT_DIR}/input/ghsl/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V2_0.tif",
    output:
        ghs_mod_resampled = "{OUTPUT_DIR}/mark_paper/high_resolution/ghsl/{DATASET}/GHS_SMOD_E2020.tif",
        rwi_resampled = "{OUTPUT_DIR}/mark_paper/high_resolution/rwi/{DATASET}/rwi.tif"
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.ghs_mod_resampled})
        mkdir --parents $(dirname {output.rwi_resampled})
        
        # pull out bounding box coords into bash array
        declare -a "COORDS=($(jq '.extracts[0].bbox | @sh' {input.json_file}))"

        # Pull reference pixel resolution
        tr=$(gdalinfo {input.pop_file} | grep "Pixel Size" | head -n1 | sed -e 's/.*(//' -e 's/).*//' | awk -F, '{{print $1, $2}}')

        # Clip and resample ghs-mod
        gdalwarp \
            -te ${{COORDS[@]}} \
            -t_srs EPSG:4326 \
            -r nearest \
            -tr $tr \
            {input.ghs_mod} \
            {output.ghs_mod_resampled}
        
        # Resample rwi
        gdalwarp \
            -t_srs EPSG:4326 \
            -r nearest \
            -tr $tr \
            {input.rwi_file} \
            {output.rwi_resampled}
        """
""" 
Test with
snakemake -c1 results/mark_paper/high_resolution/ghsl/kenya-latest/GHS_SMOD_E2020.tif
"""

rule relative_risk_flood_3as:
    """
    This rule calculates the relative risk per grid cell according to a vulnerability curve.
    The two curves that we use are either the JRC residential (average of Asia, Africa, and SAmerica)
    or the Bernhofen curve (refugee paper)
    """
    input:
        flood_file="{OUTPUT_DIR}/input/{HAZARD_SLUG}/trimmed/{DATASET}/jrc_global_flood_RP{RP}.tif"
    output:
        risk_file="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/{HAZARD_SLUG}/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP{RP}.tif"
    wildcard_constraints:
        VULN_CURVE="BERN|JRC",
        RP="10|20|50|75|100|200|500"
    script:
        "./relative_flood_risk.py"
""" 
Test with
snakemake -c1 results/mark_paper/high_resolution/relative_risk/hazard-jrc-river/kenya-latest/JRC/jrc_global_flood_RP10.tif
"""

rule average_annual_risk_jrc_3as:
    """
    This rule calculates one layer (average annual relative risk) given a set of 
    return period flood maps. This rule is written specifically for the JRC maps.
    Also calcualtes the "protected" AAR - by integrating flopros layer
    """
    input:
        flood_rp_10="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP10.tif",
        flood_rp_20="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP20.tif",
        flood_rp_50="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP50.tif",
        flood_rp_75="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP75.tif",
        flood_rp_100="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP100.tif",
        flood_rp_200="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP200.tif",
        flood_rp_500="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP500.tif",
    output:
        flood_aar="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_AAR.tif",
    wildcard_constraints:
        VULN_CURVE="BERN|JRC",
    script:
        "./average_annual_risk.py"
"""
Test with
snakemake -c1 results/mark_paper/high_resolution/relative_risk/hazard-jrc-river/kenya-latest/JRC/jrc_global_flood_AAR.tif
"""

rule rasterize_flopros_3as:
    """
    We are going to rasterize the flopros layer (using the merged protection value column).
    """
    input:
        flopros="{OUTPUT_DIR}/mark_paper/flopros/{DATASET}/flopros.shp",
        rwi_file="{OUTPUT_DIR}/mark_paper/high_resolution/rwi/{DATASET}/rwi.tif"
    output:
        rasterized_flopros="{OUTPUT_DIR}/mark_paper/high_resolution/flopros/{DATASET}/flopros.tif"
    script:
        "./rasterize_flopros.py"
"""
Test with 
snakemake -c1 results/mark_paper/high_resolution/flopros/kenya-latest/flopros.tif
""" 

rule average_annual_risk_jrc_protected_3as:
    """
    This rule calculates one layer (average annual relative risk) given a set of 
    return period flood maps. This rule is written specifically for the JRC maps.
    Also calcualtes the "protected" AAR - by integrating flopros layer
    """
    input:
        flood_rp_10="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP10.tif",
        flood_rp_20="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP20.tif",
        flood_rp_50="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP50.tif",
        flood_rp_75="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP75.tif",
        flood_rp_100="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP100.tif",
        flood_rp_200="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP200.tif",
        flood_rp_500="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_RP500.tif",
        flopros="{OUTPUT_DIR}/mark_paper/high_resolution/flopros/{DATASET}/flopros.tif",
    output:
        flood_aar="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/{VULN_CURVE}/jrc_global_flood_AAR_protected.tif",
    wildcard_constraints:
        VULN_CURVE="BERN|JRC",
    script:
        "./average_annual_risk_protected.py"
"""
Test with
snakemake -c1 results/mark_paper/high_resolution/relative_risk/hazard-jrc-river/kenya-latest/JRC/jrc_global_flood_AAR_protected.tif
"""

rule concentration_index_3as:
    """
    This rule calcualtes flood risk concentration indices for different admin regions.
    Here calculating it using the JRC damage curve...
    """
    input:
        admin_areas = "{OUTPUT_DIR}/input/admin-boundaries/{ADMIN_SLUG}.geoparquet",
        json_file="{OUTPUT_DIR}/json/{DATASET}.json",
        rwi_file="{OUTPUT_DIR}/mark_paper/high_resolution/rwi/{DATASET}/rwi.tif",
        pop_file="{OUTPUT_DIR}/mark_paper/high_resolution/ghsl/{DATASET}/GHS_POP_E2020.tif",
        urban_file="{OUTPUT_DIR}/mark_paper/high_resolution/ghsl/{DATASET}/GHS_SMOD_E2020.tif",
        risk_file="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/JRC/jrc_global_flood_AAR.tif",
        risk_protected_file="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/JRC/jrc_global_flood_AAR_protected.tif"
    output:
        regional_CI = "{OUTPUT_DIR}/mark_paper/high_resolution/concentration_index/{DATASET}/concentration_index_{ADMIN_SLUG}.gpkg",
    script:
        "./concentration_index.py"
"""
Test with
snakemake -c1 results/mark_paper/high_resolution/concentration_index/kenya-latest/concentration_index_admin-level-1.gpkg
"""

rule plot_hr_concentration_curve:
    """
    This rule plots the concentration curve for the chosen admin region.
    Here calculating it using the JRC damage curve...
    """
    input:
        admin_areas = "{OUTPUT_DIR}/input/admin-boundaries/{ADMIN_SLUG}.geoparquet",
        json_file="{OUTPUT_DIR}/json/{DATASET}.json",
        rwi_file="{OUTPUT_DIR}/mark_paper/high_resolution/rwi/{DATASET}/rwi.tif",
        pop_file="{OUTPUT_DIR}/mark_paper/high_resolution/ghsl/{DATASET}/GHS_POP_E2020.tif",
        urban_file="{OUTPUT_DIR}/mark_paper/high_resolution/ghsl/{DATASET}/GHS_SMOD_E2020.tif",
        risk_file="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/JRC/jrc_global_flood_AAR.tif",
        risk_protected_file="{OUTPUT_DIR}/mark_paper/high_resolution/relative_risk/hazard-jrc-river/{DATASET}/JRC/jrc_global_flood_AAR_protected.tif"
    output:
        figure_directory = directory("{OUTPUT_DIR}/mark_paper/high_resolution/concentration_index/{DATASET}/plots/{ADMIN_SLUG}/"),
    script:
        "./plot_concentration_curves.py"
"""
Test with
snakemake -c1 results/mark_paper/concentration_index/kenya-latest/plots/admin-level-0
"""

COUNTRIES = ["algeria-latest",
            "angola-latest",
            "burkina-faso-latest",
            "benin-latest",
            "botswana-latest",
            "burundi-latest",
            "cape-verde-latest",
            "cameroon-latest",
            "central-african-republic-latest",
            "chad-latest",
            "congo-latest",
            "comores-latest",
            "democratic-republic-congo-latest",
            "djibouti-latest",
            "egypt-latest",
            "ethiopia-latest",
            "eritrea-latest",
            "swaziland-latest",
            "equatorial-guinea-latest",
            "gabon-latest",
            "ghana-latest",
            "guinea-latest",
            "guinea-bissau-latest",
            "ivory-coast-latest",
            "kenya-latest",
            "lesotho-latest",
            "libya-latest",
            "liberia-latest",
            "madagascar-latest",
            "mali-latest",
            "malawi-latest",
            "mauritania-latest",
            "mauritius-latest",
            "morocco-latest",
            "mozambique-latest",
            "namibia-latest",
            "nigeria-latest",
            "niger-latest",
            "rwanda-latest",
            "south-africa-latest",
            "tanzania-latest",
            "tunisia-latest",
            "togo-latest",
            "sao-tome-and-principe-latest",
            "seychelles-latest",
            "sierra-leone-latest",
            "senegal-and-gambia-latest",
            "somalia-latest",
            "south-sudan-latest",
            "sudan-latest",
            "uganda-latest",
            "zambia-latest",
            "zimbabwe-latest",
]

ADMINS = ["admin-level-0"]
OUTPUT_DIRS = ["results"]

rule all_countries:
    input:
        expand("{OUTPUT_DIR}/mark_paper/high_resolution/concentration_index/{DATASET}/concentration_index_{ADMIN_SLUG}.gpkg",
               OUTPUT_DIR=OUTPUT_DIRS, DATASET=COUNTRIES, ADMIN_SLUG=ADMINS)

