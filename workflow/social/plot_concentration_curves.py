"""
This script plots the concentration curve for the administration level of interest.

"""

import logging
import json
import os

import rasterio
from rasterio.features import geometry_mask
import pandas as pd
import geopandas as gpd
import numpy as np
import shapely
import matplotlib.pyplot as plt

if __name__ == "__main__":

    try:
        admin_path: str = snakemake.input["admin_areas"]
        bbox_path: str = snakemake.input["json_file"]
        rwi_path: str = snakemake.input["rwi_file"]
        pop_path: str = snakemake.input["pop_file"]
        urban_path: str = snakemake.input["urban_file"]
        risk_path: str = snakemake.input["risk_file"]
        risk_protected_path: str = snakemake.input["risk_protected_file"]
        output_dir: str = snakemake.output["figure_directory"]
        administrative_level: int = snakemake.wildcards.ADMIN_SLUG
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

admin_level = int(administrative_level.replace("admin-level-", ""))
logging.info(f"Plotting Concentration Curves at Admin Level {admin_level}.")

if not os.path.exists(output_dir):
    logging.info("Creating Concentration Curve Figure Directory")
    os.makedirs(output_dir)

logging.info("Reading raster data.")
with rasterio.open(rwi_path) as rwi_src, rasterio.open(pop_path) as pop_src, \
    rasterio.open(urban_path) as urban_src, rasterio.open(risk_path) as risk_src, \
    rasterio.open(risk_protected_path) as risk_protected_src:
    rwi = rwi_src.read(1)
    pop = pop_src.read(1)
    urban = urban_src.read(1)
    risk = risk_src.read(1)
    risk_protected = risk_protected_src.read(1)
    affine = risk_src.transform

logging.info("Reading level {admin_level} admin boundaries")
# Read bbox bounds
with open(bbox_path, "r") as fp:
    extracts, = json.load(fp)["extracts"]
    minx, miny, maxx, maxy = extracts["bbox"]
assert minx < maxx
assert miny < maxy
bbox = shapely.geometry.box(minx, miny, maxx, maxy)

admin_areas: gpd.GeoDataFrame = gpd.read_parquet(admin_path)
area_unique_id_col = f"GID_{admin_level}"
area_name_col = f"NAME_{admin_level}"
# retain names of larger, encompassing adminstrative units
contextual_name_cols = [f"NAME_{i}" for i in range(0, admin_level)]
admin_areas = admin_areas[[area_name_col, area_unique_id_col, *contextual_name_cols, "geometry"]]
# Clip admin file to region of interest
clipped_admin_areas = admin_areas[admin_areas.intersects(bbox)]
logging.info(f"Found {len(clipped_admin_areas)} admin areas intersecting with bbox")

logging.info("Looping over admin regions and plotting concentration curves")
results = [] # List for collecting results
 # Loop over each admin region
count = 1 # for calculating progress
for idx, region in clipped_admin_areas.iterrows():
    logging.info(f"Beginning analysis of region {count} of {len(clipped_admin_areas)}.")
    count += 1
    logging.info("Clipping rasters to admin region.")
    # Get the geometry for the current admin region
    geom = region["geometry"].__geo_interface__
    # Create a mask from the geometry
    mask_array = geometry_mask([geom],
                                transform=affine,
                                invert=True,
                                out_shape=rwi.shape)
    # Use the mask to clip each raster by setting values outside the region to nan
    rwi_clip = np.where(mask_array, rwi, np.nan)
    pop_clip = np.where(mask_array, pop, np.nan)
    urban_clip = np.where(mask_array, urban, np.nan)
    risk_clip = np.where(mask_array, risk, np.nan)
    risk_protected_clip = np.where(mask_array, risk_protected, np.nan)
    logging.info("Masking and flattening data")
    rwi[rwi==-999] = np.nan # convert -999 in RWI dataset to NaN
    # Mask out areas where not all rasters are valid
    mask = (
        ~np.isnan(pop_clip) &
        ~np.isnan(rwi_clip) &
        ~np.isnan(risk_clip) &
        ~np.isnan(risk_protected_clip) &
        ~np.isnan(urban_clip)
    )
    # Flatten data
    pop_flat = pop_clip[mask]
    rwi_flat = rwi_clip[mask]
    risk_flat = risk_clip[mask]
    risk_protected_flat = risk_protected_clip[mask]
    urban_flat = urban_clip[mask]
    # Mask out zero-populatoin cells
    valid = pop_flat > 0
    pop_flat = pop_flat[valid]
    rwi_flat = rwi_flat[valid]
    risk_flat = risk_flat[valid]
    risk_protected_flat = risk_protected_flat[valid]
    urban_flat = urban_flat[valid]

    logging.info("Preparing dataframes for plotting")
    # Going to calculate concentration indices for total, urban, and rural populations (un)protected
    df = pd.DataFrame({
        'pop': pop_flat,
        'rwi': rwi_flat,
        'flood': risk_flat,
        'urban': urban_flat,
    })
    df_protected = pd.DataFrame({
        'pop': pop_flat,
        'rwi': rwi_flat,
        'flood': risk_protected_flat,
        'urban': urban_flat,
    })
    df_rural = df.copy()
    df_urban = df.copy()
    df_protected_rural = df_protected.copy()
    df_protected_urban = df_protected.copy()
    df_rural = df_rural[df_rural['urban']<=13] # GHS-MOD values 13 and under
    df_urban = df_urban[df_urban['urban']>=21] # GHS-MOD values 21 and over
    df_protected_rural = df_protected_rural[df_protected_rural['urban']<=13] # GHS-MOD values 13 and under
    df_protected_urban = df_protected_urban[df_protected_urban['urban']>=21] # GHS-MOD values 21 and over

    region_code = region[area_name_col]
    admin_output_dir = os.path.join(output_dir, region_code)
    if not os.path.exists(admin_output_dir):
        logging.info("Creating new adminstrative directory for plots")
        os.makedirs(admin_output_dir)

    logging.info("Plotting Concentration Curves")
    # Define both functions 
    def compute_concentration_curve(data):
        # Sort by RWI ascending
        data = data.sort_values(by="rwi", ascending=True).copy()
        data['cum_pop'] = data['pop'].cumsum()
        data['flood_pop'] = data['flood'] * data['pop']  # flood burden
        data['cum_flood_pop'] = data['flood_pop'].cumsum()
        total_pop = data['pop'].sum()
        total_flood_pop = data['flood_pop'].sum()
        data['frac_pop'] = data['cum_pop'] / total_pop
        data['frac_flood'] = data['cum_flood_pop'] / total_flood_pop
        # Insert 0 at the beginning for a proper starting point at (0,0)
        x_vals = np.insert(data['frac_pop'].values, 0, 0)
        y_vals = np.insert(data['frac_flood'].values, 0, 0)
        return x_vals, y_vals

    def plot_concentration_curve(x, y, title, admin_output_dir):
        plt.plot(x, y, label="Concentration Curve")
        plt.plot([0,1], [0,1], "--", color='gray', label='Equality Line')
        plt.xlabel('Cumulative Population (Wealth Rank)')
        plt.ylabel('Cumulative Flood Risk Exposure')
        plt.title(title)
        plt.legend()
        plt.savefig(os.path.join(admin_output_dir, "%s.png" % title))
        plt.close()

    x, y = compute_concentration_curve(df)
    plot_concentration_curve(x, y, f"Unprotected_Concentration_Curve_{region_code}", admin_output_dir)
    r_x, r_y = compute_concentration_curve(df_rural)
    plot_concentration_curve(r_x, r_y, f"Unprotected_Rural_Concentration_Curve_{region_code}", admin_output_dir)
    u_x, u_y = compute_concentration_curve(df_urban)
    plot_concentration_curve(u_x, u_y, f"Unprotected_Urban_Concentration_Curve_{region_code}", admin_output_dir)
    p_x, p_y = compute_concentration_curve(df_protected)
    plot_concentration_curve(p_x, p_y, f"Protected_Concentration_Curve_{region_code}", admin_output_dir)
    p_r_x, p_r_y = compute_concentration_curve(df_protected_rural)
    plot_concentration_curve(p_r_x, p_r_y, f"Protected_Rural_Concentration_Curve_{region_code}", admin_output_dir)
    p_u_x, p_u_y = compute_concentration_curve(df_protected_urban)
    plot_concentration_curve(p_u_x, p_u_y, f"Protected_Urban_Concentration_Curve_{region_code}", admin_output_dir)

logging.info("Done.")