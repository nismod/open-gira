import logging
import os

import geopandas as gpd
import numpy as np


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )

    df = gpd.read_parquet(snakemake.input.raw)

    logging.info("Relabelling track_id")
    df["track_id"] = \
        df["sample"].map(lambda x: f"S{x:03d}") \
        + df["year"].map(lambda x: f"Y{x:04d}") \
        + df["tc_number"].map(lambda x: f"N{x:03d}") \

    logging.info("Wrapping longitudes to -180, 180")
    df["lon"] = df.geometry.x
    df["lat"] = df.geometry.y
    df.lon = np.where(df.lon > 180, df.lon - 360, df.lon)
    df = gpd.GeoDataFrame(
        data=df, geometry=gpd.points_from_xy(df["lon"], df["lat"], crs=4326)
    )
    df = df.drop(columns=["lon", "lat"])

    logging.info("Writing out to disk")
    os.makedirs(os.path.dirname(snakemake.output.processed), exist_ok=True)
    df.to_parquet(snakemake.output.processed)
