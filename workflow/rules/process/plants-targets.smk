"""Process powerplants and targets for each country

"""

pyfile = os.path.join(
    WORKFLOW_DIR, "scripts", "processing", "process_power_world_countries.py"
)
# Issues: "GRL"


rule process_all:
    input:
        expand(
            os.path.join(DATA_DIR, "processed", "{code}_plants.csv"),
            code=COUNTRY_CODES,
        ),
        expand(
            os.path.join(DATA_DIR, "processed", "{code}_targets.csv"),
            code=COUNTRY_CODES,
        ),


rule process_plants_targets:
    output:
        os.path.join(DATA_DIR, "processed", "{code}_plants.csv"),
        os.path.join(DATA_DIR, "processed", "{code}_targets.csv"),
    shell:
        "python3 " + pyfile + " {wildcards.code} {output[0]} {output[1]}"
