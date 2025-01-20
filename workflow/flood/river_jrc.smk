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


rule download_river_jrc:
    output:
        zip="{OUTPUT_DIR}/input/hazard-river-jrc/raw/floodMapGL_rp{RP}y.zip"
    shell:
        """
        output_dir=$(dirname {output.zip})

        wget -q -nc \
            --directory-prefix=$output_dir \
            https://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/FLOODS/GlobalMaps/floodMapGL_rp{wildcards.RP}y.zip
        """

rule extract_river_jrc:
    input:
        tiff="{OUTPUT_DIR}/input/hazard-river-jrc/raw/floodMapGL_rp{RP}y.zip"
    output:
        tiff="{OUTPUT_DIR}/input/hazard-river-jrc/raw/floodMapGL_rp{RP}y.tif"
    shell:
        """
        output_dir=$(dirname {output.tiff})
        unzip $output_dir/floodMapGL_rp{wildcards.RP}y.zip floodMapGL_rp{wildcards.RP}y.tif -d $output_dir
        """

rule all_river_jrc:
    input:
        tiffs=expand("results/input/hazard-river-jrc/raw/floodMapGL_rp{RP}y.tif", RP=[10, 20, 50, 100, 200, 500])

"""
Current step is to run this explicitly:

    snakemake -c 10 -- all_river_jrc

Before asking for any risk results, so the checkpoint download_hazard_datasets
will pick up the TIFFs in the directory.
"""
