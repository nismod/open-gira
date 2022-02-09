# Download hazard data from the file_list.txt files in each hazard directory

rule download_hazard_datasets:
    input:
        "{OUTPUT_DIR}/input/{HAZARD_SLUG}/file_list.txt"
    output:
        "{OUTPUT_DIR}/input/{HAZARD_SLUG}/raw/{FILENAME}.tif"
    shell:
        """
        mkdir --parents {wildcards.OUTPUT_DIR}/input/{wildcards.HAZARD_SLUG}/raw
        cd {wildcards.OUTPUT_DIR}/input/{wildcards.HAZARD_SLUG}/raw &&
        wget --no-clobber -i ../file_list.txt
        """

# Not sure how to test this within snakemake.
# It can be done from the command line by requiring a specific file, e.g.:
# snakemake --cores all data/input/hazard-aqueduct-coast/raw/inuncoast_historical_nosub_hist_rp0001_5.tif
