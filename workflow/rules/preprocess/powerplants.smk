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
