"""
Rule to document countries to be excluded (no .tif file)
"""


# TODO: can we infer which countries to exclude from our other inputs?
rule process_exclude_countries:
    input:
        rules.download_gadm_levels.output.admin_bounds_global_layer_per_level,
        rules.download_population_all.input.population_raster_by_country,
    params:
        output_dir=config["output_dir"],
    output:
        os.path.join(config["output_dir"], "power_processed", "exclude_countries.txt"),
    script:
        os.path.join("..", "..", "scripts", "process", "exclude_countries.py")
