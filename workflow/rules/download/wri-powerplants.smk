"""
Download WRI powerplants database

Reference
---------
https://www.wri.org/research/global-database-power-plants
"""


rule download_powerplants:
    output:
        csv = "{OUTPUT_DIR}/input/powerplants/global_power_plant_database.csv"
    shell:
        """
        mkdir -p {wildcards.OUTPUT_DIR}/input/powerplants
        cd {wildcards.OUTPUT_DIR}/input/powerplants
        wget https://wri-dataportal-prod.s3.amazonaws.com/manual/global_power_plant_database_v_1_3.zip
        unzip -o global_power_plant_database_v_1_3.zip
        """

"""
Test with:
snakemake -c1 -- results/input/powerplants/global_power_plant_database.csv
"""