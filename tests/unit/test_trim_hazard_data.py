import os
import sys

from . import runner


def test_trim_hazard_data():
    runner.run_test(
        "trim_hazard_data",
        "snakemake results/input/hazard-aqueduct-river/djibouti-latest -j1 --keep-target-files",
    )
