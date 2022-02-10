# Download hazard data from the file_list.txt files in each hazard directory
from re import sub

# This is a checkpoint because it outputs to a directory we want to ensure is up to date in later rules.
# Specifically, trim_hazard_data.smk looks for all *.tif files in the output directory.
checkpoint download_hazard_datasets:
    output:
        directory("{OUTPUT_DIR}/input/{HAZARD_SLUG}/raw")
    run:
        input_file_key = sub("^hazard-", "", wildcards.HAZARD_SLUG)
        input_file = config['hazard_datasets'][input_file_key]
        os.system((
            f"mkdir --parents {wildcards.OUTPUT_DIR}/input/{wildcards.HAZARD_SLUG}/raw &&\\"
            f"cd {wildcards.OUTPUT_DIR}/input/{wildcards.HAZARD_SLUG}/raw &&\\"
            f"wget --no-clobber -i {input_file}"
        ))

"""
Test with:
snakemake --cores all results/input/hazard-aqueduct-coast/raw/inuncoast_historical_nosub_hist_rp0001_5.tif
"""
