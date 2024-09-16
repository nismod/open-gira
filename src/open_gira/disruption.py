import geopandas as gpd
import pandas as pd
import rasterio

import snail.intersection


def filter_edges_by_raster(
    edges: gpd.GeoDataFrame,
    raster_path: str,
    failure_threshold: float,
    band: int = 1
) -> gpd.GeoDataFrame:
    """
    Remove edges from a network that are exposed to gridded values in excess of
    a given threshold.

    Args:
        edges: Network edges to consider. Must contain geometry column containing linestrings.
        raster_path: Path to raster file on disk, openable by rasterio.
        failure_threshold: Edges experiencing a raster value in excess of this will be
            removed from the network.
        band: Band of raster to read data from.

    Returns:
        Network without edges experiencing raster values in excess of threshold.
    """

    # we will create a new edge_id column, fail if we would overwrite
    assert "edge_id" not in edges.columns

    # split out edges without geometry as snail will not handle them gracefully
    print("Parition edges by existence of geometry...")
    no_geom_edges = edges.loc[edges.geometry.type.isin({None}), :].copy()
    with_geom_edges = edges.loc[~edges.geometry.type.isin({None}), :].copy()

    # we need an id to select edges by
    with_geom_edges["edge_id"] = range(len(with_geom_edges))

    print("Read raster transform...")
    grid = snail.intersection.GridDefinition.from_raster(raster_path)

    print("Prepare linestrings...")
    edges = snail.intersection.prepare_linestrings(with_geom_edges)

    print("Split linestrings...")
    splits = snail.intersection.split_linestrings(edges, grid)

    print("Lookup raster indices...")
    splits_with_indices = snail.intersection.apply_indices(splits, grid)

    print("Read raster data into memory...")
    with rasterio.open(raster_path) as dataset:
        raster = dataset.read(band)

    print("Lookup raster value for splits...")
    values = snail.intersection.get_raster_values_for_splits(splits_with_indices, raster)
    splits_with_values = splits_with_indices.copy()
    splits_with_values["raster_value"] = values

    print(f"Filter out edges with splits experiencing values in excess of {failure_threshold} threshold...")
    failed_splits_mask = splits_with_values.raster_value > failure_threshold
    failed_edge_ids = set(splits_with_values[failed_splits_mask].edge_id.unique())
    surviving_edges_with_geom = edges.loc[~edges.edge_id.isin(failed_edge_ids), :]

    print("Done filtering edges...")
    return pd.concat([no_geom_edges, surviving_edges_with_geom]).sort_index().drop(columns=["edge_id"])