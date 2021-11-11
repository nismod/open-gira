#!/usr/bin/env python
# coding: utf-8
"""Split network edges along grid cells, join hazard values.
"""
import logging
import os
import sys

import geopandas
import pandas
import rasterio
import snail

from snail.core.intersections import split_linestring
from snail.core.intersections import get_cell_indices
from tqdm import tqdm


def main(
    network_edges_path,
    flood_data_path,
    output_path_geopq,
    output_path_pq,
):
    river = pandas.read_csv(os.path.join(flood_data_path, 'aqueduct_river.csv'))
    subset = river

    # Read metadata for a single raster
    with rasterio.open(os.path.join(flood_data_path, river.iloc[0].filename)) as dataset:
        raster_width = dataset.width
        raster_height = dataset.height
        raster_transform = list(dataset.transform)

    # Read network edges
    logging.info("Read edges")
    core_edges = geopandas.read_parquet(network_edges_path)

    # Split edges
    logging.info("Split edges")
    core_splits = []
    for edge in tqdm(core_edges.itertuples(), total=len(core_edges)):
        # split edge
        splits = split_linestring(
            edge.geometry,
            raster_width,
            raster_height,
            raster_transform,
        )
        # add to collection
        for s in splits:
            core_splits.append({
            'id': edge.id,
            'geometry': s
        })
    core_splits = geopandas.GeoDataFrame(core_splits)

    logging.info("Split %d edges into %d pieces", len(core_edges), len(core_splits))

    # Add cell indices
    logging.info("Find indices")
    def get_indices(geom): 
        x, y = get_cell_indices(
            geom,
            raster_width,
            raster_height,
            raster_transform)
        x = x % raster_width
        y = y % raster_height
        return [x, y]
    core_splits['cell_index'] = core_splits.geometry.progress_apply(get_indices)

    # Add depth values
    logging.info("Add depth values")
    for raster in tqdm(subset.itertuples(), total=len(subset)):
        associate_raster(
            core_splits,
            raster.key,
            os.path.join(flood_data_path, raster.filename))

    # Write data
    logging.info("Write data")
    core_splits.to_parquet(output_path_geopq)

    logging.info("Write data without geometry")
    pandas.DataFrame(core_splits.drop(columns=["geometry"])).to_parquet(
        output_path_pq
    )

    logging.info("Done.")


def associate_raster(df, key, fname, band_number=1):
    with rasterio.open(fname) as dataset:
        band_data = dataset.read(band_number)
        df[key] = df.cell_index.apply(lambda i: band_data[i[1], i[0]])


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    tqdm.pandas()
    try:
        network_edges_path = snakemake.input["network"]
        flood_data_path = snakemake.config["aqueduct_dir"]
        output_path_geopq = snakemake.output["geoparquet"]
        output_path_pq = snakemake.output["parquet"]
    except NameError:
        (
            network_edges_path,
            flood_data_path,
            output_path_geopq,
            output_path_pq,
        ) = sys.argv[1:]
    main(
        network_edges_path, flood_data_path, output_path_geopq, output_path_pq
    )
