"""
Download all files

Run all of the download rules, whether or not any further rules in the pipeline require the data.
"""


rule download_all:
    input:
        out_fixed,
        out_events,
        out_population,
        out_GDP,
        out_powerplant,
        out_gridfinder,
        ADMIN_BOUNDS_FILE_PER_COUNTRY,
        ADMIN_BOUNDS_GLOBAL_SINGLE_LAYER,
        ADMIN_BOUNDS_GLOBAL_LAYER_PER_LEVEL,
