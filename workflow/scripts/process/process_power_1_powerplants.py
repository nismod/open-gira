"""Download the plants data to csv files
"""
import geopandas
import pandas
import pygeos


if __name__ == "__main__":
    input_file = snakemake.input.powerplants  # type: ignore
    output_file = snakemake.output.powerplants  # type: ignore
    global_boxes_path = snakemake.input.global_boxes  # type: ignore
    box_id = snakemake.wildcards.BOX  # type: ignore

    boxes = geopandas.read_file(global_boxes_path) \
        .set_index("box_id")
    box = boxes.loc[[f"box_{box_id}"], :]
    box = box.set_crs("epsg:4326")

    powerplants = pandas.read_csv(
        input_file,
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
        columns={"gppd_idnr": "source_id"}
    )
    powerplants["type"] = "source"
    powerplants["geometry"] = pygeos.points(powerplants.longitude, powerplants.latitude)
    powerplants = geopandas.GeoDataFrame(powerplants.drop(columns=["longitude","latitude"]))
    powerplants = powerplants.set_crs("epsg:4326")
    powerplants = powerplants.sjoin(box).rename(columns={'index_right': box.index.name})
    powerplants.to_parquet(output_file)
