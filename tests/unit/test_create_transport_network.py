import os
import sys

import common

sys.path.insert(0, os.path.dirname(__file__))


def test_create_transport_network():
    common.run_test(
        'create_transport_network',
        (
            'snakemake '
            'results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_nodes.geoparquet '
            'results/geoparquet/tanzania-latest_filter-highway-core/slice-0_road_edges.geoparquet '
            '-j1 --keep-target-files'
        )
    )
