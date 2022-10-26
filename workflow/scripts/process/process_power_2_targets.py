"""This file processes the target data
"""
import sys

import geopandas as gpd
import numpy as np
import rasterio
import rasterio.features
import rasterio.mask
import xarray as xr
from pyproj import Geod
from rasterstats import gen_zonal_stats
from shapely.geometry import shape
from tqdm import tqdm


def get_target_areas(targets_file):
    geod = Geod(ellps="WGS84")
    geoms = []
    areas_km2 = []

    # Targets: Binary raster showing locations predicted to be connected to distribution grid.
    with rasterio.open(targets_file) as dataset:
        # Read the dataset's valid data mask as a ndarray.
        mask = dataset.dataset_mask()
        # Extract feature shapes and values from the array.
        for geom, val in rasterio.features.shapes(mask, transform=dataset.transform):
            if val > 0:
                feature = shape(geom)
                geoms.append(feature)
                area_m2, _ = geod.geometry_area_perimeter(feature)
                areas_km2.append(area_m2 / 1e6)

    return gpd.GeoDataFrame({
        "area_km2": areas_km2,
        "geometry": geoms
    })


def get_population(targets, population_file):
    stats = gen_zonal_stats(
        targets.geometry,
        population_file,
        stats=[],
        add_stats={"nansum": np.nansum},  # count NaN as zero for summation
        all_touched=True,  # possible overestimate, but targets grid is narrower than pop
    )
    populations = [
        d["nansum"]
        for d in tqdm(stats, desc=f"Population progress", total=len(targets))
    ]

    return populations


def get_gdp_pc(targets, gdp_file):
    gdp = xr.open_dataset(gdp_file)
    centroids = targets.geometry.centroid
    gdp_pc = gdp.sel(x=centroids.x, y=centroids.y)
    return gdp_pc


if __name__ == '__main__':
    try:
        population_file = snakemake.input.population  # type: ignore
        gdp_file = snakemake.input.gdp  # type: ignore
        targets_file = snakemake.input.targets  # type: ignore
        output_file = snakemake.output.targets  # type: ignore
    except:
        population_file = sys.argv[1]
        gdp_file = sys.argv[2]
        targets_file = sys.argv[3]
        output_file = sys.argv[4]

    targets = get_target_areas(targets_file)
    targets["type"] = "target"
    targets["population"] = get_population(targets, population_file)
    targets["gdp_pc"] = get_gdp_pc(targets, gdp_file)
    targets["gdp"] = targets.population * targets.gdp_pc

    targets.to_parquet(output_file)
