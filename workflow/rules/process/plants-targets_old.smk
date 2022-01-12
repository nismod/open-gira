"""Process powerplants and targets for each country

"""

out_world_network_with_gdp = os.path.join(
    DATA_DIR, "processed", "world_network_with_gdp.gpkg"
)


def pyfile(num):
    ext = ["_countries.py", "_network.py", "_assigngdp.py"]
    return os.path.join(
        WORKFLOW_DIR,
        "scripts",
        "processing",
        "process_power_" + str(num) + ext[num - 1],
    )


rule process_all_plants_targets:
    input:
        expand(
            os.path.join(DATA_DIR, "processed", "{code}_plants.csv"),
            code=COUNTRY_CODES,
        ),
        expand(
            os.path.join(DATA_DIR, "processed", "{code}_targets.csv"),
            code=COUNTRY_CODES,
        ),


rule process_123:
    input:
        out_world_network_with_gdp,


rule process_1:
    output:
        os.path.join(DATA_DIR, "processed", "{code}_plants.csv"),
        os.path.join(DATA_DIR, "processed", "{code}_targets.csv"),
    shell:
        "python3 " + pyfile(1) + " {wildcards.code} {output[0]} {output[1]}"


rule process_2:
    input:
        expand(
            os.path.join(DATA_DIR, "processed", "{code}_plants.csv"),
            code=COUNTRY_CODES,
        ),
        expand(
            os.path.join(DATA_DIR, "processed", "{code}_targets.csv"),
            code=COUNTRY_CODES,
        ),
    output:
        expand(
            os.path.join(DATA_DIR, "processed", "world_{attr}.gpkg"),
            attr=["edges", "network", "plants", "targets"],
        ),
    shell:
        "python3 " + pyfile(2) + ' """' + str(COUNTRY_CODES) + '"""'


rule process_3:
    input:
        expand(
            os.path.join(DATA_DIR, "processed", "world_{attr}.gpkg"),
            attr=["edges", "network", "plants", "targets"],
        ),
    output:
        out_world_network_with_gdp,
        os.path.join("data","processed","edge_gdp_sorted.txt"),
    shell:
        "python3 " + pyfile(3)
