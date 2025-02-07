"""
Given a gridded social dataset (either RWI or HDI), we will overlay with the flood raster
and return all flood exposed gridded social cells.

Flood map is resampled to the resolution of the social dataset.
"""

import logging

import rasterio
from rasterio.warp import reproject, Resampling
import numpy as np

if __name__ == "__main__":

    try:
        flood_path: str = snakemake.input["flood_file"]
        rwi_path: str = snakemake.input["rwi_file"]
        output_path: str = snakemake.output["rwi_exposure"]
        return_period: int = snakemake.wildcards["RP"]
        hazard_dataset: str = snakemake.wildcards["HAZARD_SLUG"]
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Calculating social flood exposure for {hazard_dataset} at RP {return_period}.")

logging.info("Reading raster data.")
with rasterio.open(rwi_path) as social_src:
    social = social_src.read(1)
    social_meta = social_src.meta

with rasterio.open(flood_path) as flood_src:
    flood = flood_src.read(1)
    flood_transform = flood_src.transform

    logging.info("Resampling flood data")
    resampled_array = np.empty_like(social, dtype=np.float32)

    # User rasterios reproject to resample flood raster
    reproject(
        source=flood,
        destination=resampled_array,
        src_transform=flood_src.transform,
        src_crs=flood_src.crs,
        dst_transform=social_meta["transform"],
        dst_crs=social_meta["crs"],
        resampling=Resampling.max # using max resampling here, as downsampling... could also use average
    )

logging.info("Performing flood exposure analysis")
depth_threshold = 0 # using a 0 m flood depth exposure threshold
exposed_social = np.where(resampled_array>depth_threshold, social, 0)

logging.info("Writing output raster.")
with rasterio.open(output_path, "w", **social_meta) as dst:
    dst.write(exposed_social.astype(np.float32), 1)

logging.info("Done.")