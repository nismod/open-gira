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
