#!/usr/bin/env python
# coding: utf-8
"""
Create an image showing exposed road length by combining exposure tifs,
coastline data, and administrative boundary information.
"""
import logging
import geopandas as gp
import rasterio.plot
import os
import matplotlib.pyplot as plt
import shapely.geometry as shape
import re
import glob

if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    try:
        hazard_dir = snakemake.input["hazard_dir"]
        coastline = os.path.join(
            snakemake.input["coastline"], snakemake.params["coastline"]
        )
        boundaries = os.path.join(
            snakemake.input["boundaries"], snakemake.params["boundaries"]
        )
        output_dir = snakemake.output[0]

        try:
            opts_dict = snakemake.config["exposure_tifs"]["plot"]
        except KeyError:
            opts_dict = {}

    except NameError:
        print(sys.argv)
        (hazard_dir, coastline, boundaries, output_dir, opts_dict) = sys.argv[1:]
        # hazard_dir = '../../results/exposure/tanzania-latest_filter-highway-core/hazard-aqueduct-river/'
        # coastline = '../../results/input/coastlines/ne_10m_ocean/ne_10m_ocean.shp'
        # boundaries = '../../results/input/admin-boundaries/ne_50m/ne_50m_admin_0_countries.shp'
        # output_dir = '../../results/exposure/tanzania-latest_filter-highway-core/hazard-aqueduct-river/img/'
        # opts_dict = {}

    # Load up the options from the opts_dict
    def opt(s, default=None):
        return opts_dict[s] if s in opts_dict.keys() else default

    raster_opts = opt("raster", {"cmap": "Reds"})
    coastline_opts = opt("coastline", {"facecolor": "#87cefa"})
    boundary_opts = opt("boundary", {"edgecolor": "#000000"})

    logging.info(f"Saving images to {output_dir}")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info(f"Created directory {output_dir}")

    logging.debug("Preparing coastline data.")
    coast = gp.read_file(coastline)

    logging.debug("Preparing administrative boundary data.")
    bounds = gp.read_file(boundaries)

    hazard_files = glob.glob(os.path.join(hazard_dir, "*.tif"))
    logging.info(f"Found {len(hazard_files)} .tif files in {hazard_dir}")

    for hazard_tif in hazard_files:
        with rasterio.open(hazard_tif) as hazard:
            # To plot everything in the same way and in the same place
            # we need to make sure everything uses the same
            # Coordinate Reference System (CRS).
            # The CRS we'll use as the boss is the hazard raster file CRS.
            plt_crs = hazard.crs.to_epsg()
            band = hazard.read()
            hazard_rect = shape.Polygon(
                [
                    (hazard.bounds.left, hazard.bounds.top),
                    (hazard.bounds.right, hazard.bounds.top),
                    (hazard.bounds.right, hazard.bounds.bottom),
                    (hazard.bounds.left, hazard.bounds.bottom),
                ]
            )

            # Make axes
            fig, ax = plt.subplots(dpi=300)

            # Plot raster
            logging.debug("Plotting raster data.")
            rasterio.plot.show(hazard, ax=ax, zorder=1, **raster_opts)

            # Plot coastline
            coast_tmp = coast.to_crs(plt_crs)
            coast_tmp = coast_tmp.geometry.clip(hazard_rect)
            coast_tmp.plot(ax=ax, edgecolor="none", zorder=2, **coastline_opts)

            # Plot boundaries
            bounds_tmp = bounds.to_crs(plt_crs)
            bounds_tmp = bounds_tmp.geometry.clip(hazard_rect)
            bounds_tmp.plot(ax=ax, facecolor="none", zorder=3, **boundary_opts)

            output_path = os.path.join(
                output_dir, re.sub("\\.tiff?$", ".png", os.path.basename(hazard_tif))
            )
            logging.info(f"Saving {output_path}")

            # debugging an image file comparison failure in test_make_exposure_img.py
            # TODO: remove this block
            import matplotlib
            import sys
            print(f"{matplotlib.get_backend()=}", file=sys.stderr)

            plt.savefig(output_path)
