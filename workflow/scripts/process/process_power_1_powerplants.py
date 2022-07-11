"""
This file downloads the plants data to csv files
"""

from importing_modules import *
from process_power_functions import idxbox

try:
    output_dir = snakemake.params["output_dir"]
except:
    output_dir = sys.argv[1]

if __name__ == "__main__":

    # preliminary function variables
    with open(
        os.path.join(output_dir, "power_processed", "world_boxes_metadata.txt"), "r"
    ) as filejson:
        world_boxes_metadata = json.load(filejson)
    boxlen = world_boxes_metadata["boxlen"]
    lat_max = world_boxes_metadata["lat_max"]
    lon_min = world_boxes_metadata["lon_min"]
    num_cols = world_boxes_metadata["num_cols"]
    tot_boxes = world_boxes_metadata["tot_boxes"]

    """returns all powerplants (processed data)"""
    powerplant_file = os.path.join(
        output_dir, "input", "powerplants", "global_power_plant_database.csv"
    )
    powerplants = pd.read_csv(
        powerplant_file, dtype={"other_fuel3": object}
    )  # dtype added to suppress error

    powerplants["box_id"] = idxbox(
        powerplants["latitude"],
        powerplants["longitude"],
        boxlen,
        lat_max,
        num_cols,
        lon_min,
        tot_boxes,
    )  # sort for box

    powerplants["geometry"] = powerplants.apply(
        lambda row: shape(
            {"type": "Point", "coordinates": [row.longitude, row.latitude]}
        ),
        axis=1,
    )

    cols = [
        "gppd_idnr",
        "name",
        "capacity_mw",
        "estimated_generation_gwh_2017",
        "primary_fuel",
        "box_id",
        "geometry",
    ]
    powerplants = gpd.GeoDataFrame(powerplants[cols])

    powerplants.rename(columns={"gppd_idnr": "source_id"}, inplace=True)

    powerplants["type"] = "source"

    for box_id, powerplants_box in tqdm(
        powerplants.groupby("box_id"),
        desc="saving powerplant csv",
        total=len(powerplants["box_id"].unique()),
    ):
        all_boxes_path = os.path.join(
            output_dir, "power_processed", "all_boxes", f"{box_id}"
        )
        if not os.path.exists(all_boxes_path):
            os.makedirs(all_boxes_path)
        p = os.path.join(all_boxes_path, f"powerplants_{box_id}.csv")
        powerplants_box.to_csv(p, index=False)

    with open(
        os.path.join(output_dir, "power_processed", "world_boxes_metadata.txt"), "r"
    ) as filejson:
        tot_boxes = json.load(filejson)["tot_boxes"]

    cols[0] = "source_id"
    cols[-1] = "type"
    cols.append("geometry")
    for id in tqdm(
        range(int(float(tot_boxes))), desc="empty folders", total=float(tot_boxes)
    ):  # create empty ones
        all_boxes_pp_file = os.path.join(
            output_dir,
            "power_processed",
            "all_boxes",
            f"box_{id}",
            f"powerplants_box_{id}.csv",
        )
        if not os.path.exists(all_boxes_pp_file):
            gpd.GeoDataFrame(columns=cols).to_csv(all_boxes_pp_file, index=False)
