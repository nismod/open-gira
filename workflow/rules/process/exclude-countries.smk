"""Rule to document countries to be excuded (no .tif file)

"""


rule process_excludecountries:
    input:
        out_adminboundaries,
        out_population,
    output:
        os.path.join("data", "adminboundaries", "exclude_countries.txt"),
    shell:
        "python3 " + os.path.join(
        "workflow", "scripts", "process", "exclude_countries.py"
        )
