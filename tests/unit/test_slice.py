import os
import sys
import common

sys.path.insert(0, os.path.dirname(__file__))


def test_slice():
    common.run_test(
        "slice",
        "snakemake results/slices/tanzania-mini_filter-highway-core -j1 --keep-target-files",
    )
