import logging
import os

import geopandas as gpd
import shapely


def subset_file_by_intersection(geoms: gpd.GeoDataFrame, in_path: str, out_path: str) -> None:
    """
    Subset some geoparquet data on disk by a set of geometries and write to
    disk again.

    Args:
        geoms: Geometries of interest
        in_path: Location of dataset to subset
        out_path: Location to write subset
    """

    logging.info(f"Subsetting {in_path}")
    df: gpd.GeoDataFrame = gpd.read_parquet(in_path)
    subset: gpd.GeoDataFrame = df.sjoin(geoms, how="inner")

    logging.info(f"Writing subset to {out_path}")
    subset[df.columns].to_parquet(out_path)
    return


if __name__ == "__main__":

    gridfinder_path: str = snakemake.input.gridfinder
    targets_path: str = snakemake.input.targets
    powerplants_path: str = snakemake.input.powerplants
    admin_bounds_path: str = snakemake.input.admin_bounds
    gridfinder_out_path: str = snakemake.output.gridfinder
    targets_out_path: str = snakemake.output.targets
    powerplants_out_path: str = snakemake.output.powerplants
    country_iso_a3: str = snakemake.wildcards.COUNTRY_ISO_A3

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    # read admin bounds for country in question
    countries = gpd.read_parquet(admin_bounds_path).rename(columns={"GID_0": "iso_a3"})
    country = countries[countries.iso_a3 == country_iso_a3]
    country = country[["iso_a3", "geometry"]].copy()

    os.makedirs(os.path.dirname(gridfinder_out_path), exist_ok=True)

    subset_file_by_intersection(country, gridfinder_path, gridfinder_out_path)
    subset_file_by_intersection(country, targets_path, targets_out_path)
    subset_file_by_intersection(country, powerplants_path, powerplants_out_path)
