"""
Given a flood depth map and a specified vulnerability curve, this script calculates the relative risk (0-1)
for each floodded grid cell.
"""

import logging

import rasterio
import numpy as np

if __name__ == "__main__":

    try:
        flood_path: str = snakemake.input["flood_file"]
        output_path: str = snakemake.output["risk_file"]
        return_period: int = snakemake.wildcards["RP"]
        hazard_dataset: str = snakemake.wildcards["HAZARD_SLUG"]
        vuln_dataset: str = snakemake.wildcards["VULN_CURVE"]
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Calculating flood risk for {hazard_dataset} at RP {return_period} using {vuln_dataset} vulnerability curve.")

logging.info(f"Initializing vulnerability curve")
vuln_curves = {"JRC": ([0, 0.5, 1, 1.5, 2, 3, 4, 5, 6], [0, 0.35, 0.53, 0.66, 0.77, 0.89, 0.94, 0.98, 1]), # average of Asia, sAmerica, Africa Residential curve
             "BERN": ([0, 0.15, 0.5, 1.5], [0, 0.33, 0.67, 1])
             }
vuln_curve = vuln_curves[vuln_dataset]

logging.info("Reading raster data.")

with rasterio.open(flood_path) as flood_src:
    flood = flood_src.read(1)
    flood_transform = flood_src.transform
    flood_meta = flood_src.meta

logging.info("Performing flood risk analysis")

def vectorized_damage(depth, heights, damage_percents):
    '''
    Vectorized damage function
    Apply damage function given a flood depth.
    Function also needs as input the damage function heights > damage_percents
    '''
    # Use np.interp for vectorized linear interpolation
    damage_percentage = np.interp(depth, heights, damage_percents)
    return damage_percentage

risk = vectorized_damage(flood, vuln_curve[0], vuln_curve[1])

logging.info("Writing output raster.")
with rasterio.open(output_path, "w", **flood_meta) as dst:
    dst.write(risk.astype(np.float32), 1)

logging.info("Done.")