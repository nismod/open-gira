import os
import sys

from . import runner


def test_slice():
    runner.run_test(
        "slice",
        "snakemake results/slices/djibouti-latest_filter-road -j1 --keep-target-files",
    )
