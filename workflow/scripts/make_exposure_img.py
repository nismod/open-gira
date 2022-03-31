#!/usr/bin/env python
# coding: utf-8
"""
Create an image showing exposed road length by combining exposure tifs,
coastline data, and administrative boundary information.
"""
import logging
import geopandas as gp
import rasterio
import rasterio.plot
import cartopy.crs as crs
import os
import matplotlib.pyplot as plt
import shapely.geometry as shape


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    # move to settings later
    epsg = 3857
    try:
        hazard_tif = snakemake.input['tif'][0]
        coastline = snakemake.input['coastline']
        boundaries = snakemake.input['boundaries']
        output_path = snakemake.output[0]
        opts_dict = snakemake.config['exposure_tifs']
    except NameError:
        # print(sys.argv)
        # (
        #     hazard_tif,
        #     coastline,
        #     boundaries,
        #     output_path,
        #     opts_dict
        # ) = sys.argv[1:]
        hazard_tif = '../../results/exposure/tanzania-mini/hazard-aqueduct-river/exposure_inunriver_rcp4p5_MIROC-ESM-CHEM_2030_rp00100.tif'
        coastline = '../../results/input/coastlines/ne_10m_ocean/ne_10m_ocean.shp'
        boundaries = '../../results/input/admin-boundaries/ne_50m/ne_50m_admin_0_countries.shp'
        output_path = '../../results/exposure/tanzania-mini/hazard-aqueduct-river/img/exposure_inunriver_rcp4p5_MIROC-ESM-CHEM_2030_rp00100.png'
        opts_dict = {}

    # Load up the options from the opts_dict

    logging.info(f"Generating image {output_path}")
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
        logging.info(f"Created directory {os.path.dirname(output_path)}")

    with rasterio.open(hazard_tif) as hazard:
        # To plot everything in the same way and in the same place
        # we need to make sure everything uses the same
        # Coordinate Reference System (CRS).
        # The CRS we'll use as the boss is the hazard raster file CRS.
        plt_crs = hazard.crs.to_epsg()
        band = hazard.read()
        hazard_rect = shape.Polygon([
            (hazard.bounds.left, hazard.bounds.top),
            (hazard.bounds.right, hazard.bounds.top),
            (hazard.bounds.right, hazard.bounds.bottom),
            (hazard.bounds.left, hazard.bounds.bottom),
        ])

        # Make axes
        fig, ax = plt.subplots(dpi=300)

        # Plot raster
        logging.debug("Plotting raster data.")
        rasterio.plot.show(hazard, ax=ax, zorder=3)

        # Plot coastline
        logging.debug("Plotting coastline data.")
        coast = gp.read_file(coastline).to_crs(plt_crs)
        coast = coast.geometry.clip(hazard_rect)
        coast.plot(ax=ax, edgecolor='none', facecolor='#87cefa', zorder=1)

        # Plot boundaries
        logging.debug("Plotting administrative boundary data.")
        bounds = gp.read_file(boundaries).to_crs(plt_crs)
        bounds = bounds.geometry.clip(hazard_rect)
        bounds.plot(ax=ax, edgecolor='#333333', facecolor='none', zorder=2)

        plt.savefig(output_path)
