import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_join_data():
    common.run_test(
        'join_data',
        (
            'snakemake results/tanzania-mini_filter-highway-core_hazard-aqueduct-river.geoparquet '
            '-j1 --keep-target-files'
        )
    )
