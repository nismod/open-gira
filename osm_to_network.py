#!/usr/bin/env python
# coding: utf-8
"""Read OSM pbf, load driving network, write to parquet as-is and core only.
"""
import logging
import os
import sys

import pyrosm


def main(pbf_path, outputs_path):
    osm = pyrosm.OSM(pbf_path)
    slug = os.path.basename(pbf_path).replace(".osm.pbf", "")
    nodes, edges = osm.get_network(nodes=True, network_type="driving")

    # Write direct from pyrosm driving
    logging.info("Nodes: %d", len(nodes))
    nodes.to_parquet(
        os.path.join(
            outputs_path,
            f'{slug}-roads-edges.geoparquet'))
    logging.info("Edges: %d", len(edges))
    edges.to_parquet(
        os.path.join(
            outputs_path,
            f'{slug}-roads-nodes.geoparquet'))

    core = (
        'motorway_link',
        'motorway',
        'trunk_link',
        'trunk',
        'primary_link',
        'primary',
        'secondary_link',
        'secondary',
        'tertiary_link',
        'tertiary',
    )
    core_edges = edges[edges.highway.isin(core)]
    logging.info("Core edges: %d", len(core_edges))

    select_columns = [
        'bridge', 'highway', 'lanes', 'maxspeed', 'oneway',
        'smoothness', 'surface', 'tracktype', 'tunnel', 'width',
        'id', 'name', 'osm_type', 'geometry', 'u', 'v', 'length'
    ]
    core_edges = core_edges[select_columns]
    core_edges.to_parquet(
        os.path.join(
            outputs_path,
            f'{slug}-roads_core-edges.geoparquet'))
    logging.info("Done.")


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

    # Ignore geopandas parquet implementation warnings
    # NB though that .geoparquet is not the format to use for archiving.
    import warnings
    warnings.filterwarnings('ignore', message='.*initial implementation of Parquet.*')

    logging.info("Start")
    pbf_path, outputs_path = sys.argv[1:]
    main(pbf_path, outputs_path)
