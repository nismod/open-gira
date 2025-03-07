"""
This script runs the concentration index analysis at the given administrative level.

Available administrative layers in GADM:
    - National (level 0)
    - State/province/equivalent (level 1)
    - County/district/equivalent (level 2)
    - Smaller Level 3 or 4.

"""

import logging
import json

import rasterio
from rasterio.features import geometry_mask
import pandas as pd
import geopandas as gpd
import numpy as np
import shapely

if __name__ == "__main__":

    try:
        admin_path: str = snakemake.input["admin_areas"]
        bbox_path: str = snakemake.input["json_file"]
        rwi_path: str = snakemake.input["rwi_file"]
        pop_path: str = snakemake.input["pop_file"]
        urban_path: str = snakemake.input["urban_file"]
        risk_path: str = snakemake.input["risk_file"]
        risk_protected_path: str = snakemake.input["risk_protected_file"]
        output_path: str = snakemake.output["regional_CI"]
        administrative_level: int = snakemake.wildcards.ADMIN_SLUG
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

admin_level = int(administrative_level.replace("admin-level-", ""))
logging.info(f"Calculating concentration indices at Admin Level {admin_level}.")

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

logging.info("Looping over admin regions and calculating concentration indices")
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

    logging.info("Preparing dataframes for analysis")
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

    logging.info("Calculating Concentration Index")
    # Define function
    def calculate_CI(df):
        # Sort dataframe by wealth
        df = df.sort_values(by="rwi", ascending=True)
        # Calculate cumulative population rank (to represent distribution of people)
        df['cum_pop'] = df['pop'].cumsum()
        # Calculate total pop of sample
        total_pop = df['pop'].sum()
        if total_pop == 0:
            return np.nan
        # Calculate fractional rank of each row
        df['rank'] = (df['cum_pop'] - 0.5*df['pop']) / total_pop
        try:
            # Calcualte weighted mean of flood risk
            weighted_mean_flood = np.average(df['flood'], weights=df['pop'])
        except ZeroDivisionError:
            return np.nan
        if weighted_mean_flood == 0:
            return np.nan
        # Calculate weighted sum of (flood * rank * pop)
        sum_xR = (df['flood'] * df['rank'] * df['pop']).sum()
        # Calculate Concentration Index
        CI = (2 * sum_xR) / (df['pop'].sum() * weighted_mean_flood) - 1
        return CI

    CI = calculate_CI(df)
    CI_protected = calculate_CI(df_protected)
    CI_urban = calculate_CI(df_urban)
    CI_rural = calculate_CI(df_rural)
    CI_protected_urban = calculate_CI(df_protected_urban)
    CI_protected_rural = calculate_CI(df_protected_rural)

    results.append({
        area_unique_id_col: region[area_unique_id_col],
        area_name_col: region[area_name_col],
        "CI": CI,
        "CI_protected": CI_protected,
        "CI_urban": CI_urban,
        "CI_rural": CI_rural,
        "CI_protected_urban": CI_protected_urban,
        "CI_protected_rural": CI_protected_rural,
        "Population": np.sum(pop_flat),
        "rwi_count": np.count_nonzero(rwi_flat),
        "geometry": region["geometry"]
    })

logging.info("Writing reults to GeoPackage.")
results_gdf = gpd.GeoDataFrame(results, geometry="geometry")
results_gdf.to_file(output_path, driver="GPKG")

logging.info("Done.")