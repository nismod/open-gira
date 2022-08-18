import os
import sys

from . import runner


def test_download_hazard_datasets():
    runner.run_test(
        "download_hazard_datasets",
        "snakemake results/input/hazard-aqueduct-river/raw -j1 --keep-target-files",
    )
