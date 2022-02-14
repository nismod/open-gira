import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_download_hazard_datasets():
    common.run_test(
        'download_hazard_datasets',
        'snakemake results/input/hazard-aqueduct-river/raw -j1 --keep-target-files'
    )
