"""
Given a gridded social dataset (either RWI or HDI), we will overlay with the population dataset
and calculate the population weighted social metric at the chosen administrative region.

Code is based loosely on the following repo https://github.com/worldbank/RWI 

Available administrative layers in GADM:
    - National (level 0)
    - State/province/equivalent (level 1)
    - County/district/equivalent (level 2)
    - Smaller Level 3 or 4.

Seldom does a country have boundaries for every level.
"""

import logging
import json

import rasterio
import geopandas as gpd
from rasterstats import zonal_stats
import numpy as np
import shapely

WGS84_EPSG = 4326

if __name__ == "__main__":

    try:
        admin_areas_path: str = snakemake.input["admin_areas"]
        social_grid_path: str = snakemake.input["rwi_file"]
        bbox_path: str = snakemake.input["bbox_file"]
        pop_grid_path: str = snakemake.input["pop_file"]
        admin_level_slug: int = snakemake.wildcards.ADMIN_SLUG
        pop_weight_path: int = snakemake.output["pop_weighted_rwi"]
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info("Reading raster data.")
with rasterio.open(social_grid_path) as social_src, rasterio.open(pop_grid_path) as pop_src:
    pop = pop_src.read(1)
    social = social_src.read(1)
    # Get the affine transformation matrix
    affine = pop_src.transform

logging.info("Multiplying rasters")
# Debug issue with rwi no data value
social[social == -999] = np.nan
raster_product = pop * social


admin_level = int(admin_level_slug.replace("admin-level-", ""))
logging.info("Reading level {admin_level} admin boundaries.")

# Read bbox bounds
with open(bbox_path, "r") as fp:
    extracts, = json.load(fp)["extracts"]
    minx, miny, maxx, maxy = extracts["bbox"]
assert minx < maxx
assert miny < maxy
bbox = shapely.geometry.box(minx, miny, maxx, maxy)

admin_areas: gpd.GeoDataFrame = gpd.read_parquet(admin_areas_path)
area_unique_id_col = f"GID_{admin_level}"
area_name_col = f"NAME_{admin_level}"
# retain names of larger, encompassing adminstrative units
contextual_name_cols = [f"NAME_{i}" for i in range(0, admin_level)]
admin_areas = admin_areas[[area_name_col, area_unique_id_col, *contextual_name_cols, "geometry"]]
# Clip admin file to region of interest
clipped_admin_areas = admin_areas[admin_areas.intersects(bbox)]
logging.info(f"Found {len(clipped_admin_areas)} admin areas intersecting with bbox")

logging.info("Running zonal stats.")
# Ensure CRS compatibility 
raster_product_crs = pop_src.crs.to_dict()
clipped_admin_areas = clipped_admin_areas.to_crs(raster_product_crs)
# Run zonal stats on raster product and population totals
admin_raster_product = zonal_stats(clipped_admin_areas, raster_product, affine=affine, stats=["sum"], nodata=-999)
admin_population = zonal_stats(clipped_admin_areas, pop, affine=affine, stats=["sum"], nodata=-999)

logging.info("Looping through admin regions to calculate regional population weights.")
# Initialize list for storing population weighted social metric
regional_pop_weighted_social_list = []
# Initialize list for storing population 
regional_pop_list = []

# Loop through each region
for region_raster_product, region_population in zip(admin_raster_product, admin_population):
    # Total product (pop * social) in region
    region_product_sum = region_raster_product["sum"]
    # Total population in region
    region_pop_sum = region_population["sum"]

    # Debug: ensure None is treated as 0
    if region_product_sum is None:
        region_product_sum = 0
    if region_pop_sum is None:
        region_pop_sum = 0

    # Calculate population weighted social metric for the region
    if region_pop_sum != 0:
        population_weighted_social = region_product_sum / region_pop_sum
        regional_pop_weighted_social_list.append(population_weighted_social)
    else:
        # If total population is zero, assign NaN
        regional_pop_weighted_social_list.append(np.nan)

    # Append total population to list
    regional_pop_list.append(region_pop_sum)

logging.info("Storing results to dataframe")
clipped_admin_areas['Population_Weighted_RWI'] = regional_pop_weighted_social_list
clipped_admin_areas['Total_Population'] = regional_pop_list

logging.info("Writing social population weighted admin file to disk.")
clipped_admin_areas.to_parquet(pop_weight_path)

logging.info("Done")

