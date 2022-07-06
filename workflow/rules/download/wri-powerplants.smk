"""Download WRI powerplants database

Reference
---------
https://www.wri.org/research/global-database-power-plants
"""

out_powerplant = os.path.join(
    config["output_dir"], "input", "powerplants", "global_power_plant_database.csv"
)


rule download_powerplants:
    output:
        out_powerplant,
    shell:
        f"""
        mkdir -p {config['output_dir']}/input/powerplants
        cd {config['output_dir']}/input/powerplants
        wget https://wri-dataportal-prod.s3.amazonaws.com/manual/global_power_plant_database_v_1_3.zip
        unzip -o global_power_plant_database_v_1_3.zip
        """
