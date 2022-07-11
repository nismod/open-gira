import os
import sys

import common

sys.path.insert(0, os.path.dirname(__file__))


def test_join_network_nodes():
    common.run_test(
        'join_network_nodes',
        (
            'snakemake results/tanzania-latest_filter-highway-core/road_nodes.geoparquet '
            '-j1 --keep-target-files'
        )
    )
