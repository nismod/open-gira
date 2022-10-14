"""Download the plants data to csv files
"""
import json
import os
import sys
import warnings

import geopandas as gpd
import pandas as pd
import pygeos
from tqdm import tqdm

from process_power_functions import idxbox

warnings.filterwarnings("ignore", category=DeprecationWarning)


if __name__ == "__main__":
    try:
        input_file = snakemake.input.powerplants  # type: ignore
        output_file = snakemake.output.powerplants  # type: ignore
    except:
        input_file = sys.argv[1]
        output_file = sys.argv[2]

    powerplants = pd.read_csv(
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
    powerplants = gpd.GeoDataFrame(powerplants.drop(columns=["longitude","latitude"]))
    powerplants.to_parquet(output_file)
