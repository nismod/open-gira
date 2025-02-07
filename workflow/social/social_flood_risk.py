"""
This script extracts the RWI grid cells exposed to flooding at various depth-vulnerability thresholds.
It returns raster datasets of the gridded exposed RWI data.
Raster analysis performed at the resolution of the RWI dataset (2.4km). As such, the flood data is 
resampled to this resolution. 
Using the flood vulnerability thresholds from Bernhofen et al. https://iopscience.iop.org/article/10.1088/1748-9326/acd8d0/meta 
"""

import logging

import rasterio
from rasterio.warp import reproject, Resampling
import numpy as np
import pandas as pd

if __name__ == "__main__":

    try:
        flood_path: str = snakemake.input["flood_file"]
        rwi_path: str = snakemake.input["rwi_file"]
        return_period: int = snakemake.wildcards["RP"]
        hazard_dataset: str = snakemake.wildcards["HAZARD_SLUG"]
        vhigh_output: str = snakemake.output["rwi_vhigh_risk"]
        high_output: str = snakemake.output["rwi_high_risk"]
        medium_output: str = snakemake.output["rwi_medium_risk"]
        low_output: str = snakemake.output["rwi_low_risk"]
        exposure_output: str = snakemake.output["rwi_exposure"]
    except NameError:
        raise ValueError("Must be run via snakemake.")

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Calculating social flood risk for {hazard_dataset} at RP {return_period}.")

logging.info("Initializing vulnerability thresholds")
"""
These thresholds are from Bernhofen et al. (2023) and relate flood depth to "level" of risk.
0 m = no risk; 0-0.15 m = low risk; 0.15-0.5 m = medium risk; 0.5-1.5 m = high risk; >1.5 = very high risk
"""
vulnerability_thresholds = {
    "very_high":(1.5, 100),
    "high": (0.5, 1.5),
    "medium": (0.15, 0.5),
    "low": (0, 0.15)
}

logging.info("Reading raster data.")
with rasterio.open(rwi_path) as rwi_src:
    rwi = rwi_src.read(1)
    rwi_meta = rwi_src.meta

with rasterio.open(flood_path) as flood_src:
    flood = flood_src.read(1)
    flood_transform = flood_src.transform
    flood = np.where(flood>0, flood, 0) # remove NaN values

    logging.info("Resampling flood data")
    resampled_array = np.empty_like(rwi, dtype=np.float64)

    # User rasterio's reproject to resample flood raster
    reproject(
        source=flood,
        destination=resampled_array,
        src_transform=flood_src.transform,
        src_crs=flood_src.crs,
        dst_transform=rwi_meta["transform"],
        dst_crs=rwi_meta["crs"],
        resampling=Resampling.average # using average resampling here
    )

    logging.info("Performing flood risk analysis")
    for risk_level in vulnerability_thresholds:
        logging.info(f"Working on risk level: {risk_level}.")
        min_depth = vulnerability_thresholds[risk_level][0]
        max_depth = vulnerability_thresholds[risk_level][1]

        exposed_social = np.where((resampled_array>min_depth) & (resampled_array<=max_depth), rwi, np.nan)

        logging.info("Writing output risk raster")
        # This is quite long (TODO: come up with a better way so we can automatically assign output paths)
        if risk_level == "very_high":
            with rasterio.open(vhigh_output, "w", **rwi_meta) as dst:
                dst.write(exposed_social.astype(np.float32), 1)
        if risk_level == "high":
            with rasterio.open(high_output, "w", **rwi_meta) as dst:
                dst.write(exposed_social.astype(np.float32), 1)
        if risk_level == "medium":
            with rasterio.open(medium_output, "w", **rwi_meta) as dst:
                dst.write(exposed_social.astype(np.float32), 1)
        if risk_level == "low":
            with rasterio.open(low_output, "w", **rwi_meta) as dst:
                dst.write(exposed_social.astype(np.float32), 1)

    logging.info("Performing flood exposure analysis")
    exposed_social = np.where(resampled_array>0, rwi, np.nan)

    logging.info("Writing output exposure raster")
    with rasterio.open(exposure_output, "w", **rwi_meta) as dst:
        dst.write(exposed_social.astype(np.float64), 1)

logging.info("Done.")
