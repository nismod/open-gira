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

from pyproj import Geod
from snail.core.intersections import get_cell_indices, split_linestring
from tqdm import tqdm


def main(network_edges_path, attrs, hazard_data_path, hazard_data_csv, outputs_path):
    # Filename to use for output
    network_slug = os.path.basename(network_edges_path).replace(".geoparquet", "")
    hazard_slug = os.path.basename(hazard_data_csv).replace(".csv", "")
    slug = f"{network_slug}_{hazard_slug}"

    # Read hazard metadata
    # This is a config/steering file for this script, assumes hazards are all on the same
    # grid, and hazards and networks are in the same CRS.
    hazards = pandas.read_csv(hazard_data_csv)

    # Read metadata for a single raster
    with rasterio.open(os.path.join(hazard_data_path, hazards.iloc[0].filename)) as dataset:
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
            split_data = {
                'geometry': s
            }
            for attr in attrs:
                split_data[attr] = getattr(edge, attr)
            core_splits.append(split_data)
    core_splits = geopandas.GeoDataFrame(core_splits)

    logging.info("Split %d edges into %d pieces", len(core_edges), len(core_splits))

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

    logging.info("Segment length")
    geod = Geod(ellps="WGS84")
    core_splits['length_km'] = core_splits.geometry.progress_apply(geod.geometry_length) / 1e3

    logging.info("Add hazard values")
    for raster in tqdm(hazards.itertuples(), total=len(hazards)):
        associate_raster(
            core_splits,
            raster.key,
            os.path.join(hazard_data_path, raster.filename))

    logging.info("Write data")
    core_splits.to_parquet(os.path.join(outputs_path, f'{slug}.splits.geoparquet'))

    logging.info("Write data without geometry")
    pandas.DataFrame(core_splits.drop(columns=['geometry'])) \
        .to_parquet(os.path.join(outputs_path, f'{slug}.splits.parquet'))

    logging.info("Done.")


def associate_raster(df, key, fname, band_number=1):
    with rasterio.open(fname) as dataset:
        band_data = dataset.read(band_number)
        df[key] = df.cell_index.apply(lambda i: band_data[i[1], i[0]])


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    tqdm.pandas()
    print(sys.argv)
    network_edges_path, attrs, hazard_data_path, hazard_csv, outputs_path = sys.argv[1:]
    attrs = attrs.split(",")
    main(network_edges_path, attrs, hazard_data_path, hazard_csv, outputs_path)
