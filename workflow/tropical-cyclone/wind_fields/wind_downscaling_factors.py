"""
Estimate downscaling factors from gradient to surface level winds, given surface
roughness raster.
"""

import logging
import os

import numpy as np
import rasterio
import rioxarray
import xarray as xr

from open_gira.wind import power_law_scale_factors
from open_gira.wind_plotting import plot_downscale_factors

if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    # wind speed altitudes

    # gradient level clearly not realistic, but we vary it to fit our estimated wind
    # speeds to observations (or to better model results)
    # the 18m figure is as a result of minimising pixel-wise errors between this model
    # and that used in Done et al. 2020 with a physical boundary layer
    GRADIENT_LEVEL_METRES = 18
    SURFACE_LEVEL_METRES = 10

    logging.info("Calculating downscaling factors from surface roughness raster")

    # surface roughness raster for downscaling winds with
    try:
        surface_roughness_raster: xr.DataArray = rioxarray.open_rasterio(snakemake.input.surface_roughness)
    except rasterio.errors.RasterioIOError:
        logging.info("Surface roughness raster is empty, writing empty downscaling factors matrix.")

        # assure (empty) files exist
        np.save(snakemake.output.downscale_factors, np.array([]))
        os.system(f"touch {snakemake.output.downscale_factors_plot}")

        sys.exit(0)

    # (1, y, x) where 1 is the number of surface roughness bands
    # select the only value in the band dimension
    surface_roughness: np.ndarray = surface_roughness_raster[0].values

    # calculate factors to scale wind speeds from gradient-level to surface level,
    # taking into account surface roughness as defined by the raster
    downscaling_factors = power_law_scale_factors(
        surface_roughness,
        SURFACE_LEVEL_METRES,
        GRADIENT_LEVEL_METRES
    )

    logging.info("Saving downscaling factors to disk")
    np.save(snakemake.output.downscale_factors, downscaling_factors)

    logging.info("Mapping downscaling factors and saving to disk")
    plot_downscale_factors(
        downscaling_factors,
        "Wind downscaling factors",
        snakemake.output.downscale_factors_plot
    )

    logging.info("Done estimating downscaling factors")