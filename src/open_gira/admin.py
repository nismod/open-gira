"""
Tools for manipulating administrative region data.
"""

import logging

import pandas as pd
import geopandas as gpd
import shapely


def merge_gadm_admin_levels(preference: pd.DataFrame, alternative: pd.DataFrame) -> pd.DataFrame:
    """
    Geospatial data is often aggregated at admin 'levels', as exemplified by the
    GADM project. These integer levels are 0 for nation states, 1 for the
    regions of the state, 2 for sub-regions and so on. This function is designed
    to take two different levels, e.g. 2 and 1 as `preference` and
    `alternative`, respectively and combine them. Where a country has
    representation in `preference`, these rows will be returned by default. If
    there is a country in `alternative` that is not in `preference`, this
    country's rows from `alternative` will be merged into the returned result.

    Args:
        preference: Table with rows to return by default. Must contain an ISO_A3
            column.
        alternative: Table with rows to choose from if `preference` is missing a
            country. Must contain an ISO_A3 column.

    Returns:
        Table of rows representing the union of countries in the passed tables,
        with data from `preference` by default, but `alternative` if `preference`
        data is not available.
    """

    substitute_countries = set(alternative.ISO_A3) - set(preference.ISO_A3)
    logging.info(f"Gap filling with: {substitute_countries}")

    substitute_regions: pd.DataFrame = alternative[alternative["ISO_A3"].isin(substitute_countries)]

    merged = pd.concat([preference, substitute_regions])

    return merged.sort_values("ISO_A3")


def get_administrative_data(file_path: str, to_epsg: int = None) -> gpd.GeoDataFrame:
    """
    Read administrative data (country ISO, country geometry) from disk

    Arguments:
        file_path (str): Location of file with country data:
            containing an ISO three letter code as 'GID_0' and a geometry as
            'geometry'
        to_epsg (int): EPSG code to project data to

    Returns:
        gpd.GeoDataFrame: Table of country and geometry data with:
            'iso_a3' and 'geometry' columns
    """

    # read file
    gdf = gpd.read_file(file_path)

    # check schema is as expected
    expected_columns = {"GID_0", "geometry"}
    assert expected_columns.issubset(set(gdf.columns.values))

    # reproject if desired
    if to_epsg is not None:
        gdf = gdf.to_crs(epsg=to_epsg)

    # rename these columns first so we don't have to do this twice (to nodes and edges) later
    gdf.rename(columns={"GID_0": "iso_a3"}, inplace=True)

    # subset, sort by iso_a3 and return
    return gdf[["iso_a3", "geometry"]].sort_values(by=["iso_a3"], ascending=True)


def boundary_geom(gdf: gpd.GeoDataFrame, iso_a3: str) -> shapely.Geometry:
    """Given administrative data, return the boundary geometry for a given ISO3 country code
    """
    return gdf.set_index("iso_a3").loc[iso_a3, "geometry"]
