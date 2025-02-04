"""Download and extract JRC river flood return period maps

https://data.jrc.ec.europa.eu/dataset/jrc-floods-floodmapgl_rp50y-tif

The global river flood hazard maps are a gridded data set representing
inundation along the river network, for seven different flood return periods
(from 1-in-10-years to 1-in-500-years). The input river flow data for the new
maps are produced by means of the open-source hydrological model LISFLOOD, while
inundation simulations are performed with the hydrodynamic model LISFLOOD-FP.
The extent comprises the entire world with the exception of Greenland and
Antarctica and small islands with river basins smaller than 500km2.

Cell values indicate water depth (in m). The maps can be used to assess the
exposure of population and economic assets to river floods, and to perform flood
risk assessments. The dataset is created as part of the Copernicus Emergency
Management Service. NOTE: this dataset is not an official flood hazard map (for
details and limitations please refer to related publications).

Citation:

Baugh, Calum; Colonese, Juan; D'Angelo, Claudia; Dottori, Francesco; Neal,
Jeffrey; Prudhomme, Christel; Salamon, Peter (2024): Global river flood hazard
maps. European Commission, Joint Research Centre (JRC) [Dataset] PID:
http://data.europa.eu/89h/jrc-floods-floodmapgl_rp50y-tif
"""

rule download_jrc_flood:
    output:
        zip="{OUTPUT_DIR}/input/jrc_flood/floodMapGL_rp{RP}y.zip"
    shell:
        """
        output_dir=$(dirname {output.zip})

        wget -q -nc \
            --directory-prefix=$output_dir \
            https://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/FLOODS/GlobalMaps/floodMapGL_rp{wildcards.RP}y.zip
        """
rule extract_jrc_flood:
    input:
        tiff="{OUTPUT_DIR}/input/jrc_flood/floodMapGL_rp{RP}y.zip"
    output:
        tiff="{OUTPUT_DIR}/input/jrc_flood/floodMapGL_rp{RP}y.tif"
    shell:
        """
        output_dir=$(dirname {output.tiff})
        unzip $output_dir/floodMapGL_rp{wildcards.RP}y.zip floodMapGL_rp{wildcards.RP}y.tif -d $output_dir
        """

rule all_jrc_flood:
    input:
        tiffs=expand("results/input/jrc_flood/floodMapGL_rp{RP}y.tif", RP=[10, 20, 50, 100, 200, 500])

rule merge_jrc_flood:
    """
    Merges all the jrc flood tiles
    """
    input:
        raw_folder="{OUTPUT_DIR}/input/{HAZARD_SLUG}/raw/"
    output:
        RP10="{OUTPUT_DIR}/input/{HAZARD_SLUG}/merged/jrc_global_flood_RP10.tif",
        RP20="{OUTPUT_DIR}/input/{HAZARD_SLUG}/merged/jrc_global_flood_RP20.tif",
        RP50="{OUTPUT_DIR}/input/{HAZARD_SLUG}/merged/jrc_global_flood_RP50.tif",
        RP75="{OUTPUT_DIR}/input/{HAZARD_SLUG}/merged/jrc_global_flood_RP75.tif",
        RP100="{OUTPUT_DIR}/input/{HAZARD_SLUG}/merged/jrc_global_flood_RP100.tif",
        RP200="{OUTPUT_DIR}/input/{HAZARD_SLUG}/merged/jrc_global_flood_RP200.tif",
        RP500="{OUTPUT_DIR}/input/{HAZARD_SLUG}/merged/jrc_global_flood_RP500.tif"
    script:
        "./merge_jrc.py"
"""
Test with
snakemake -c1 results/input/hazard-jrc-river/merged/jrc_global_flood_RP10.tif
"""
    


