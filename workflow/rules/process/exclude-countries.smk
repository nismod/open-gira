"""Rule to document countries to be excluded (no .tif file)

"""


rule process_exclude_countries:
    input:
        out_adminboundaries_levels,
        out_population,
    output:
        os.path.join("data", "adminboundaries", "exclude_countries.txt"),
    script:
        os.path.join("..", "..", "scripts", "process", "exclude_countries.py")
