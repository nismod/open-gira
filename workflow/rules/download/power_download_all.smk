"""
Download all files

Run all of the download rules, whether or not any further rules in the pipeline require the data.
"""


rule download_all:
    input:
        STORMS_RETURN_PERIOD,
        STORMS_EVENTS,
        rules.download_population_all.input.population_raster_by_country,
        rules.download_GDP.output.gdp_datasets,
        rules.download_powerplants.output.powerplants_global,
        rules.download_gridfinder.output.electricity_grid_global,
        rules.download_gadm.output.admin_bounds_global_single_layer,
        rules.download_gadm_levels.output.admin_bounds_global_layer_per_level,
        rules.download_gadm_all_countries.input.admin_bounds_file_per_country,
