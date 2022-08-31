"""Rule to document countries to be excluded (no .tif file)

"""


rule process_exclude_countries:
    input:
        ADMIN_BOUNDS_GLOBAL_LAYER_PER_LEVEL,
        out_population,
    params:
        output_dir=config["output_dir"],
    output:
        os.path.join(config["output_dir"], "power_processed", "exclude_countries.txt"),
    script:
        os.path.join("..", "..", "scripts", "process", "exclude_countries.py")
