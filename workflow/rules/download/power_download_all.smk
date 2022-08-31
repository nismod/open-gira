"""
Download all files

Run all of the download rules, whether or not any further rules in the pipeline require the data.
"""


rule download_all:
    input:
        STORMS_RETURN_PERIOD,
        STORMS_EVENTS,
        POPULATION_RASTER_BY_COUNTRY,
        GDP_DATASETS,
        POWERPLANTS_GLOBAL,
        ELECTRICITY_GRID_GLOBAL,
        ADMIN_BOUNDS_FILE_PER_COUNTRY,
        ADMIN_BOUNDS_GLOBAL_SINGLE_LAYER,
        ADMIN_BOUNDS_GLOBAL_LAYER_PER_LEVEL,
