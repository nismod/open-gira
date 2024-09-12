"""
Plotting utilities.
"""


import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.ops import split


def chop_at_antimeridian(gdf: gpd.GeoDataFrame, drop_null_geometry: bool = False) -> gpd.GeoDataFrame:
    """
    Cut LineStrings either side of antimeridian, then drop the fragments that
        intersect with antimeridian.

    Warning: Will create new rows (split geometries) with duplicate indices.

    Args:
        gdf: Table with geometry to chop at antimeridian
        drop_null_geometry: If true, drop any null geometry rows before plotting

    Returns:
        Table, potentially with new rows. No rows in the table should have
            geometries that cross the antimeridian.
    """
    if drop_null_geometry:
        gdf = gdf.loc[~gdf.geometry.isna(), :]

    assert set(gdf.geometry.type) == {'LineString'}

    def split_on_meridian(gdf: gpd.GeoDataFrame, meridian: shapely.geometry.LineString) -> gpd.GeoDataFrame:
        return gdf.assign(geometry=gdf.apply(lambda row: split(row.geometry, meridian), axis=1)).explode(index_parts=False)

    xlim = 179.9
    ylim = 90

    split_e = split_on_meridian(gdf, shapely.geometry.LineString([(xlim, ylim), (xlim, -ylim)]))
    split_e_and_w = split_on_meridian(split_e, shapely.geometry.LineString([(-xlim, ylim), (-xlim, -ylim)]))

    def crosses_antimeridian(row: pd.Series) -> bool:
        """
        Check if there are longitudes in a geometry that are near the antimeridian
            (i.e. -180) and both sides of it. If so, return true.
        """
        x, _ = row.geometry.coords.xy
        longitudes_near_antimeridian = np.array(x)[np.argwhere(np.abs(np.abs(x) - 180) < xlim).ravel()]
        if len(longitudes_near_antimeridian) == 0:
            return False
        hemispheres = np.unique(np.sign(longitudes_near_antimeridian))
        if (-1 in hemispheres) and (1 in hemispheres):
            return True
        else:
            return False

    return split_e_and_w[~split_e_and_w.apply(crosses_antimeridian, axis=1)]


def figure_size(
    min_x: float,
    min_y: float,
    max_x: float,
    max_y: float,
    max_plot_width_in: float = 16,
    max_plot_height_in: float = 9
) -> tuple[float, float]:
    """
    Given bounding box, calculate size of figure in inches for equal aspect ratio.
    Will not exceed `max_plot_width_in` or `max_plot_height_in`.

    Example usage:
    When used in conjunction with a GeoDataFrame named `gdf`:
    ```
    f, ax = plt.subplots(figsize=figure_size(*gdf.total_bounds))
    ```

    Args:
        min_x: Minimum x value in data coordinates
        min_y: Minimum y value in data coordinates
        max_x: Maximum x value in data coordinates
        max_y: Maximum y value in data coordinates

    Returns:
        width in inches, height in inches
    """
    x_span = max_x - min_x
    y_span = max_y - min_y
    aspect_ratio = y_span / x_span
    max_plot_width_in = 16
    max_plot_height_in = 9

    if max_plot_width_in * aspect_ratio < max_plot_height_in:
        # tall
        x_in = max_plot_width_in
        y_in = max_plot_width_in * aspect_ratio

    else:
        # wide
        x_in = max_plot_height_in / aspect_ratio
        y_in = max_plot_height_in

    return x_in, y_in
