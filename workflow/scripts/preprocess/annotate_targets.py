"""
Annotate electricity consuming areas with encompassed GDP and population.
"""

import logging

import geopandas as gpd
import pandas as pd
import xarray as xr


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
    admin_path: str = snakemake.input.admin
    population_path: str = snakemake.input.population
    gdp_path: str = snakemake.input.gdp
    targets_path: str = snakemake.input.targets
    output_path: str = snakemake.output.targets

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    logging.info("Reading targets file")
    targets = gpd.read_parquet(targets_path)
    logging.info(f"Read {len(targets)} targets from disk")

    logging.info("Annotating population per target")
    population = pd.read_csv(population_path)
    targets = targets.merge(population, on="id").rename(columns={"population_sum": "population"})

    logging.info("Extracting GDP per target")
    # ~3 minutes CPU time for the globe
    targets["gdp_pc"] = annotate_gdp_pc(targets, gdp_path)

    # if we can't find a GDP per capita figure, set it to zero
    # we require a number for every target, so we can allocate power from sources to sinks
    nan_gdp_mask = targets["gdp_pc"].isna()
    logging.info(f"Setting GDP per capita to 0 for {len(targets[nan_gdp_mask])} targets")
    targets.loc[nan_gdp_mask, "gdp_pc"] = 0

    targets["gdp"] = targets.population * targets.gdp_pc

    logging.info("Tagging targets with country ISO alpha-3 code")
    admin = gpd.read_parquet(admin_path, columns=["GID_0", "geometry"])

    # a naive targets.sjoin(admin) using polygons and intersects predicate will take ~7 hours
    # instead, use point-in-polygon, i.e. target representative point _within_ admin polygon
    # should now take a few seconds

    # create copy of targets with representative points rather than full polygon geometry
    rep_targets = targets.copy()
    rep_targets.geometry = targets.geometry.representative_point()

    # do the join to find the containing country (or at least, one of them for border cases)
    target_country_join = rep_targets.sjoin(admin, predicate="within")
    target_country_join = target_country_join[["id", "GID_0"]].rename(columns={"GID_0": "iso_a3"})

    # merge the iso_a3 column back into the polygon geodataframe
    targets = targets.merge(target_country_join, on="id")

    targets["asset_type"] = "target"

    logging.info("Writing targets to disk")
    targets.to_parquet(output_path)
