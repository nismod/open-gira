"""
Annotate electricity consuming areas with encompassed GDP and population.
"""

import logging
import multiprocessing
import os

import numpy as np
import geopandas as gpd
import pandas as pd
import rasterio
import rasterio.features
import rasterio.mask
import xarray as xr
from rasterstats import gen_zonal_stats


def annotate_population(targets: gpd.GeoDataFrame, population_path: str) -> gpd.GeoDataFrame:
    """
    Given a `targets` table with polygonal geometries, lookup population
    encompassed by each from `population_path`.

    Args:
        targets: Table containing geometries to get population for
        population_path: Path to population raster file
    """

    with rasterio.open(population_path) as dataset:
        crs = dataset.crs.data
    stats = gen_zonal_stats(
        targets.to_crs(crs).geometry,  # reprojected for raster
        population_path,
        stats=[],
        add_stats={"nansum": np.nansum},  # count NaN as zero for summation
        all_touched=True,  # possible overestimate, but targets grid is narrower than pop
    )
    # fill masked values, in case of nodata
    populations = np.ma.array([d["nansum"] for d in stats]).filled(fill_value=0)
    populations = np.nan_to_num(populations.astype(float))

    targets["population"] = populations
    return targets


def annotate_gdp_pc(targets: gpd.GeoDataFrame, gdp_path: str) -> pd.Series:
    """
    Lookup GDP per capita for each target geometry from a netCDF file.

    Args:
        targets: Table of target geometries
        gdp_path: Path to netCDF file containing `GDP_per_capita_PPP` with
            `longitude` and `latitude` coordinates.

    Returns:
        Series of GDP per capita figures
    """

    ds = xr.open_dataset(gdp_path)

    centroids = targets.geometry.centroid
    df = pd.DataFrame({"x": centroids.x, "y": centroids.y})

    df["gdp_pc"] = ds.GDP_per_capita_PPP.sel(
        longitude=df.x.to_xarray(),
        latitude=df.y.to_xarray(),
        time=2015,
        method="nearest",
    )
    return df.gdp_pc


if __name__ == '__main__':
    population_path: str = snakemake.input.population
    gdp_path: str = snakemake.input.gdp
    targets_path: str = snakemake.input.targets
    output_path: str = snakemake.output.targets
    parallel: bool = snakemake.config["process_parallelism"]

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    logging.info("Reading targets file")
    targets = gpd.read_parquet(targets_path)
    targets["asset_type"] = "target"

    logging.info("Extracting population per target")
    # N.B. requires ~200 minutes of CPU time for globe
    targets["id"] = range(len(targets))
    if parallel:
        chunked_pop = []
        chunk_size: int = np.ceil(len(targets) / os.cpu_count()).astype(int)
        args = [(targets.iloc[i: i + chunk_size, :].copy(), population_path) for i in range(0, len(targets), chunk_size)]
        with multiprocessing.Pool() as pool:
            chunked_pop = pool.starmap(annotate_population, args)
        targets = pd.concat(chunked_pop).sort_values("id")
    else:
        targets = annotate_population(targets, population_path)

    logging.info("Extracting GDP per target")
    # ~3 minutes CPU time for the globe
    targets["gdp_pc"] = annotate_gdp_pc(targets, gdp_path)
    targets["gdp"] = targets.population * targets.gdp_pc

    logging.info("Writing targets to disk")
    targets.to_parquet(output_path)
