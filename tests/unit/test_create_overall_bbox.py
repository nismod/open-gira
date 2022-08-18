import os
import sys

from . import runner


def test_create_overall_bbox():
    runner.run_test(
        "create_overall_bbox",
        "snakemake results/json/djibouti-latest.json -j1 --keep-target-files",
    )
