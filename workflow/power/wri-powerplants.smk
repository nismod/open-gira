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

rule parse_powerplants:
    """
    Parse powerplant data for world and save in convenient format
    """
    input:
        powerplants=rules.download_powerplants.output.csv
    output:
        powerplants="{OUTPUT_DIR}/power/powerplants.geoparquet",
    run:
        import geopandas as gpd
        import pandas as pd

        powerplants = pd.read_csv(
            input.powerplants,
            # dtype for other_fuel3 added to suppress error
            dtype={"other_fuel3": object},
            usecols=(
                "gppd_idnr",
                "name",
                "capacity_mw",
                "estimated_generation_gwh_2017",
                "primary_fuel",
                "longitude",
                "latitude",
            )
        ).rename(
            columns={
                "gppd_idnr": "source_id",
                "capacity_mw": "power_mw"
            }
        )
        powerplants["asset_type"] = "source"
        powerplants["geometry"] = gpd.points_from_xy(
            powerplants.longitude,
            powerplants.latitude,
            crs="epsg:4326"
        )
        powerplants = gpd.GeoDataFrame(
            powerplants.drop(
                columns=["longitude","latitude"]
            )
        )
        powerplants.to_parquet(output.powerplants)

"""
Test with:
snakemake -c1 results/power/powerplants.geoparquet
"""
