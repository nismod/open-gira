"""
Transform land cover raster into a surface roughness raster based on wind grid
specification
"""

import logging
import os

import pandas as pd
import rioxarray
from rasterio.errors import RasterioIOError
import numpy as np

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

if __name__ == "__main__":

    try:
        land_cover = rioxarray.open_rasterio(snakemake.input.land_cover)
    except RasterioIOError:
        logging.info("Found empty land cover map, creating empty surface roughness raster...")
        os.system(f"touch {snakemake.output.surface_roughness}")
        sys.exit(0)
    wind_grid = rioxarray.open_rasterio(snakemake.input.wind_grid)
    cover_roughness = pd.read_csv(snakemake.input.land_cover_roughness_mapping, comment="#")

    # build a lookup array where the category is the index
    # and the value for a given index/category is the surface roughness in metres
    logging.info("Build land cover -> surface roughness mapping...")
    lookup_length = max(set(cover_roughness.glob_cover_2009_id)) + 1
    roughness_lookup: np.ndarray = np.zeros(lookup_length)
    for row in cover_roughness.itertuples():
        roughness_lookup[row.glob_cover_2009_id] = row.roughness_length_m

    # create surface roughness field (same shape as land cover)
    logging.info("Create surface roughness field...")
    surface_roughness_values: np.ndarray = roughness_lookup[land_cover.values]
    surface_roughness = land_cover.copy()
    surface_roughness.values = surface_roughness_values

    # surface roughness on wind grid
    logging.info("Resample to wind grid...")
    downsampled_roughness = surface_roughness.interp(x=wind_grid.x, y=wind_grid.y)

    # write out surface roughness values on wind grid
    logging.info("Save to disk...")
    downsampled_roughness.rio.to_raster(snakemake.output.surface_roughness)
