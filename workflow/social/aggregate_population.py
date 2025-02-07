"""
Given a population dataset, and the RWI dataset, this script aggregates the 
population dataset to the resolution of the RWI dataset, using an area-weighted
aggregation approach
"""

import logging

import rasterio
import geopandas as gpd
import numpy as np
import rasterstats as rs
from rasterio.features import shapes
from rasterio.transform import rowcol
from shapely.geometry import box

if __name__ == "__main__":

    try:
        pop_path: str = snakemake.input["pop_file"]
        rwi_path: str = snakemake.input["rwi_file"]
        agg_path: str = snakemake.output["agg_pop_file"]
    except NameError:
        raise ValueError("Must be run via snakemake.")

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info("Loading population raster.")
with rasterio.open(pop_path) as pop_src:
    pop = pop_src.read(1)
    pop_transform = pop_src.transform
    pop_crs = pop_src.crs
    pop_profile = pop_src.profile
    pop_nodata = pop_src.nodata

logging.info("Loading RWI rastser.")
with rasterio.open(rwi_path) as rwi_src:
    rwi_transform = rwi_src.transform
    rwi_crs = rwi_src.crs
    rwi_shape = (rwi_src.height, rwi_src.width)
    rwi_bounds = rwi_src.bounds

logging.info("Generating reference RWI grid for aggregation.")
# Get raster bounds directly from rasterio
xmin, ymin, xmax, ymax = rwi_bounds
x_coords = np.arange(xmin, xmax, rwi_transform.a)
y_coords = np.arange(ymax, ymin, -abs(rwi_transform.e))
grid_cells = []
for x in x_coords:
    for y in y_coords:
        grid_cells.append(box(x, y, x + rwi_transform.a, y + rwi_transform.e))
grid_gdf = gpd.GeoDataFrame({"geometry": grid_cells}, crs=rwi_crs)

logging.info("Performing area-weighted aggregation.")
# Note on the below. all_touched=False preserves pop totals but leads to striping
zonal_stats = rs.zonal_stats(grid_gdf, pop_path, stats="sum", all_touched=False, geojson_out=True, nodata=pop_nodata)
pop_aggregated = gpd.GeoDataFrame.from_features(zonal_stats)

logging.info("Preparing data for raster conversion")
pop_aggregated["grid_id"] = np.arange(len(pop_aggregated))
pop_raster_array = np.full(rwi_shape, 0.0)
for geom, value in zip(grid_gdf.geometry, pop_aggregated["sum"]):
    if np.isfinite(value):
        centroid = geom.centroid  # Get polygon center
        row, col = rowcol(rwi_transform, centroid.x, centroid.y)
        # Ensure the index is within raster bounds
        if 0 <= row < rwi_shape[0] and 0 <= col < rwi_shape[1]:
            pop_raster_array[row, col] = value

logging.info("Write new aggregated population raster")
pop_profile.update(dtype=rasterio.float64, transform=rwi_transform, width=rwi_shape[1], height=rwi_shape[0])
with rasterio.open(agg_path, "w", **pop_profile) as dst:
    dst.write(pop_raster_array.astype(rasterio.float64), 1)

logging.info("Done.")