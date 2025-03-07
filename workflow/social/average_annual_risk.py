"""
Given relative risk maps for flood return periods - calculate the average annual relative risk
for each floodded grid cell.
"""

import logging

import rasterio
import numpy as np

if __name__ == "__main__":

    try:
        RP10_path: str = snakemake.input["flood_rp_10"]
        RP20_path: str = snakemake.input["flood_rp_20"]
        RP50_path: str = snakemake.input["flood_rp_50"]
        RP75_path: str = snakemake.input["flood_rp_75"]
        RP100_path: str = snakemake.input["flood_rp_100"]
        RP200_path: str = snakemake.input["flood_rp_200"]
        RP500_path: str = snakemake.input["flood_rp_500"]
        aar_output_path: str = snakemake.output["flood_aar"]
        vuln_dataset: str = snakemake.wildcards["VULN_CURVE"]
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Calculating average annual relative risk using {vuln_dataset} vulnerability curve.")

logging.info("Reading raster data.")
raster_paths = [RP10_path, RP20_path, RP50_path, RP75_path, RP100_path, RP200_path, RP500_path]
flood_maps = [] # going to load rasters into this list
for path in raster_paths:
    with rasterio.open(path) as src:
        flood_maps.append(src.read(1))
        transform = src.transform
        meta = src.meta

# Create an empty array for bankful (2-year) flood
RP2 = np.zeros_like(flood_maps[0])
flood_maps.insert(0, RP2) # insert this array at beginning of list.

logging.info("Calculating annual average risk - unprotected")
RPs = np.array([2, 10, 20, 50, 75, 100, 200, 500]) # define return peridos
aep = 1 / RPs # convert to annual exceedance probability
flood_maps = np.array(flood_maps) # convert to numpy array for calculation
aar = np.trapz(flood_maps[::-1], x=aep[::-1], axis=0) # debug: need to reverse list to avoid negative AAR values

logging.info("Writing output rasters.")
with rasterio.open(aar_output_path, "w", **meta) as dst:
    dst.write(aar.astype(np.float32), 1)

logging.info("Done.")