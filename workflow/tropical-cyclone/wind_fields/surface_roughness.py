"""
Transform land cover raster into a surface roughness raster based on wind grid
specification
"""

import logging
import os
import sys

import pandas as pd
import rioxarray
from rasterio.errors import RasterioIOError
import numpy as np


logging.basicConfig(
    format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
)

if __name__ == "__main__":
    try:
        land_cover = rioxarray.open_rasterio(snakemake.input.land_cover)  # noqa: F821
    except RasterioIOError:
        logging.info(
            "Found empty land cover map, creating empty surface roughness raster..."
        )
        os.system(f"touch {snakemake.output.surface_roughness}")  # noqa: F821
        sys.exit(0)
    wind_grid = rioxarray.open_rasterio(snakemake.input.wind_grid)  # noqa: F821
    cover_roughness = pd.read_csv(
        snakemake.input.land_cover_roughness_mapping,  # noqa: F821
        comment="#",
    )

    logging.info("Build land cover -> surface roughness mapping...")
    lookup_length = max(set(cover_roughness.glob_cover_2009_id)) + 1
    roughness_lookup: np.ndarray = np.zeros(lookup_length)
    for row in cover_roughness.itertuples():
        roughness_lookup[row.glob_cover_2009_id] = row.roughness_length_m

    logging.info("Create surface roughness field...")
    input_roughness = land_cover.copy()
    input_roughness.values = roughness_lookup[land_cover.values.astype(int)]

    logging.info("Aggregate roughness to wind grid...")
    output_roughness = input_roughness.rio.reproject_match(
        wind_grid, resampling="average"
    )

    logging.info("Save to disk...")
    output_roughness.rio.to_raster(snakemake.output.surface_roughness)  # noqa: F821
