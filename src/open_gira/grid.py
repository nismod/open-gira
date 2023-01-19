import logging

import geopandas as gpd
import pandas as pd
from pyproj import Geod
import rasterio
import rasterio.features
import rasterio.mask
import shapely


def polygonise_targets(targets_path: str, extent: shapely.geometry.Polygon) -> gpd.GeoDataFrame:
    """
    Take a raster of electricity consuming 'targets' and a geometry extent and
    return a set of target polygons with computed areas.

    Args:
        targets_path: Path to raster file containing targets
        extent: Shape to mask targets by

    Returns:
        Table of target geometries and areas
    """

    geod = Geod(ellps="WGS84")
    geoms = []
    areas_km2 = []

    # Targets: Binary raster showing locations predicted to be connected to distribution grid.
    with rasterio.open(targets_path) as dataset:
        crs = dataset.crs.data

        # Read the dataset's valid data mask as a ndarray.
        try:
            box_dataset, box_transform = rasterio.mask.mask(dataset, [extent], crop=True)

            # Extract feature shapes and values from the array.
            for geom, val in rasterio.features.shapes(box_dataset, transform=box_transform):
                if val > 0:
                    feature = shapely.geometry.shape(geom)
                    geoms.append(feature)
                    area_m2, _ = geod.geometry_area_perimeter(feature)
                    areas_km2.append(abs(area_m2 / 1e6))
        except ValueError as ex:
            # could be that extent does not overlap dataset
            logging.info("Extent may not overlap targets", ex)
            pass

    return gpd.GeoDataFrame(
        data={
            "area_km2": areas_km2,
            "geometry": geoms
        },
        crs=crs
    )


def weighted_allocation(
    nodes: pd.DataFrame,
    *,
    variable_col: str,
    weight_col: str,
    component_col: str,
    asset_col: str,
    source_name: str,
    sink_name: str,
) -> pd.DataFrame:
    """
    For each network component, allocate total variable capacity of sources
    to sinks (consuming nodes), weighted by some feature of the sinks. The sum
    of the variable for all sources and sinks in a component should equal zero.

    Args:
        nodes: Should contain, `variable_col`, `weight_col`, `asset_col` and `component_col`
        variable_col: Name of column in `nodes` to distribute from `source` to `sink`
        weight_col: Name of column in `nodes` to weight allocation to sinks by
        component_col: Name of column in `nodes` to group sources and sinks by
        asset_col: Name of column in `nodes` to differentiate sources and sinks by
        source_name: Categorical for sources in `nodes[asset_col]`
        sink_name: Categorical for sinks in `nodes[asset_col]`

    Returns:
        Sink nodes with variable allocated from source nodes by weight
    """

    # find the sum of variable for each component
    c_variable_sum = nodes.loc[
        nodes[asset_col] == source_name,
        [variable_col, component_col]
    ].groupby(component_col).sum().reset_index()

    # subset to sinks
    sinks = nodes[nodes[asset_col] == sink_name]

    # find the sum of weights for each component
    c_weight_sum = sinks.loc[:, [weight_col, component_col]].groupby(component_col).sum().reset_index()

    # merge in the component sums for variable and weight
    c_variable_col = f"_component_{variable_col}"
    sinks = sinks.merge(
        c_variable_sum.rename(columns={variable_col: c_variable_col}),
        how="left",
        on=component_col
    )
    c_weight_col = f"_component_{weight_col}"
    sinks = sinks.merge(
        c_weight_sum.rename(columns={weight_col: c_weight_col}),
        how="left",
        on=component_col
    )

    # ensure every sink has a numeric (non-NaN) entry for the component sum of variable
    sinks[c_variable_col] = sinks[c_variable_col].fillna(0)

    # reallocate variable to sinks, by weight within components
    sinks[variable_col] = -1 * sinks[c_variable_col] * sinks[weight_col] / sinks[c_weight_col]

    return sinks
