"""
Routines for binning data.
"""

import geopandas as gpd
import numpy as np
import shapely


def grid_point_data(
    point_data: gpd.GeoDataFrame, col_name: str, agg_func_name: str, delta: float
) -> gpd.GeoDataFrame:
    """
    Generate regular rectilinear grid and bin GeoDataFrame point data,
    applying some aggregation function `agg_func_name` per grid cell.

    Grid will be constructed with height and width `delta`.
    Grid will snap to integer number of `delta` from CRS reference, minimally
    encompassing all `point_data`.

    Args:
        point_data: Table of data, containing `col_name` and geometry columns.
        col_name: Data column to rebin
        agg_func: Name of function to apply to multiple points within a bin, e.g. max, mean
        delta: Width and height of output raster pixels in decimal degrees
    """

    def n_cells(start: float, end: float, width: float) -> int:
        return round((end - start) / width)

    # create grid

    # unpack and snap extrema to integer number of cells from CRS reference
    min_x, min_y, max_x, max_y = map(
        lambda x: np.ceil(x / delta) * delta, point_data.total_bounds
    )

    cells = []
    for x0 in np.linspace(min_x, max_x, n_cells(min_x, max_x, delta) + 1):
        x1 = x0 + delta
        for y0 in np.linspace(min_y, max_y, n_cells(min_y, max_y, delta) + 1):
            y1 = y0 + delta
            cells.append(shapely.geometry.box(x0, y0, x1, y1))

    grid = gpd.GeoDataFrame({"geometry": cells})

    # associate points with grid
    raster = grid.sjoin(point_data, how="inner")

    # take aggregate of arbitrary number of rows (points) grouped by cell
    agg_func = getattr(raster.groupby(raster.index), agg_func_name)
    aggregated = agg_func(col_name)

    # merge geometry back in, how="right" for every grid cell, empty or not
    rebinned = aggregated.merge(
        grid[["geometry"]], left_index=True, right_index=True, how="right"
    )
    rebinned = rebinned.drop(columns=["index_right"])

    # cast to GeoDataFrame again (merge reverted to DataFrame) and return
    return gpd.GeoDataFrame(rebinned)
