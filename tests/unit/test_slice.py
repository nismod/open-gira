import os
import sys

from . import runner


def test_slice():
    runner.run_snakemake_test(
        "slice",
        (
            "results/slices/djibouti-latest_filter-road",
        )
    )
