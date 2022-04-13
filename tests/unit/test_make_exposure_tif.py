import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_make_exposure_tif():
    common.run_test(
        'make_exposure_tif',
        (
            'snakemake results/exposure/tanzania-mini_filter-highway-core/hazard-aqueduct-river/ '
            '-j1'
        )
    )
