import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_trim_hazard_data():
    common.run_test(
        'trim_hazard_data',
        'snakemake results/input/hazard-aqueduct-river/tanzania-mini -j1 --keep-target-files'
    )
