import os
import sys

from . import runner


def test_create_bbox_extracts():
    runner.run_test(
        "create_bbox_extracts",
        "snakemake results/json/djibouti-latest_extracts.geojson -j1 --keep-target-files",
    )
